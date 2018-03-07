clc;
clear;
close all;
fontSize = 20;

MAX_ITERATIONS = 3;

% ---- BALL RELATED CONSTANTS ---- % 

BALL_SIZE = 13; % pixels of the ball, CALIBRATE THE CAMERA, this comes automatically 
BALL_AREA = pi * (BALL_SIZE/2)^2;
SMALLESTACCEPTABLEAREA = BALL_SIZE*2; 
RADIUSRANGE = [round(BALL_SIZE/2-0.3*BALL_SIZE/2), round(BALL_SIZE/2+0.3*BALL_SIZE/2)]; 
ROI_WIDTH = BALL_SIZE * 8; %this can vary with the video FPS 
ROI_HEIGHT = BALL_SIZE * 4;
ADAPTIVE_WIDTH = ROI_WIDTH;
ADAPTIVE_HEIGHT = ROI_HEIGHT;

% kernel for blurring. VALUE: 3-sigma rule
BLUR_KERNEL_LENGTH = (BALL_SIZE-1)/6; 
CLOSING_KERNEL_LENGTH = round(BALL_SIZE/2); % try this value a little more

% ---- END OF BALL RELATED CONSTANTS ---- %



% ---- BEGINNING OF THE CODE ---- %

%prompt = 'Input the name of the video file: ';
%videoName = input(prompt, 's');
videoName = 'Video/orange_ball_high.mp4';
v = VideoReader(videoName);

prompt = 'Frames to be skipped: ';
skips = str2double(input(prompt, 's'));
if isnan(skips) || fix(skips) ~= skips
  disp('Please enter an integer');
  return;
end

prompt = 'Color of the ball ([o] - Orange [w] - White):  ';
ballColor = input(prompt, 's');

% --- set the initial color thresholds --- % 
if strcmp(ballColor, 'o')
    h = 45; first_threshold = 20; second_threshold = 35;
    FIRST_PASS_HSV_MIN = [(h - first_threshold)/360, 60/100, 20/100];
    FIRST_PASS_HSV_MAX = [(h + first_threshold)/360, 1, 1];
    % may be add a variable in the saturation as well 
    SECOND_PASS_HSV_MIN = [(h - second_threshold)/360, 40/100, 0];
    SECOND_PASS_HSV_MAX = [(h + second_threshold)/360, 1, 1];
elseif strcmp(ballColor, 'w')
    h = 10; first_threshold = 5; second_threshold = 12;
    FIRST_PASS_HSV_MIN = [0, 0, 0.75];
    FIRST_PASS_HSV_MAX = [1, h+first_threshold/100, 1];
    SECOND_PASS_HSV_MIN = [0, 0, 0.6];
    SECOND_PASS_HSV_MAX = [1, h+second_threshold/100, 1];
else
    disp('Invalid! Accepted Colors: [o] - Orange or [w] - White');
    return;
end


% --- var initialization --- %
skipped = 0;
positions = [];
ball_found_prev = false; 

% for accuracy (mainly low res frames) the user should select the first ROI
first_frame = true;


% --- MAIN LOOP --- %

while hasFrame(v)
    
    % skip the bad frames at the beginning
    if skipped < skips  
        if hasFrame(v)
            frame = readFrame(v);
        end
        skipped=skipped+1;
    else
        
        frame = readFrame(v);
        ball_found = 0;
        roiRect = [0, 0, 0, 0];
        candidateRoiRect = [];
        candidateBallsRoi = [];
        
        % calculate the ROI 
        if ball_found_prev == 1
            roiRect = calcRoiSize(positions(length(positions(:,1)), :), [size(frame,2), size(frame,1)], ROI_WIDTH, ROI_HEIGHT);
            roi = imcrop(frame, roiRect);
        else
            if first_frame  == true
                figure;
                imshow(frame);
                [pos(1), pos(2)] = getpts;
                roiRect = calcRoiSize(pos, [size(frame,2), size(frame,1)], ROI_WIDTH, ROI_HEIGHT);
                roi = imcrop(frame, roiRect);
                first_frame = false;
                
            else
               % IMPLEMENT EXTRAPOLATION!!
               ADAPTIVE_WIDTH = 1.4 * ADAPTIVE_WIDTH;
               ADAPTIVE_HEIGHT = 1.4 * ADAPTIVE_HEIGHT;
               % take care if not even the first ball was found 
               if isempty(positions)
                   roiRect = calcRoiSize(pos, [size(frame,2), size(frame,1)], ADAPTIVE_WIDTH, ADAPTIVE_HEIGHT);
               else
                   roiRect = calcRoiSize(positions(length(positions(:,1)), :), [size(frame,2), size(frame,1)], ADAPTIVE_WIDTH, ADAPTIVE_HEIGHT);
               end
               roi = imcrop(frame, roiRect);
            end
        end
        
        roiBlurred = imgaussfilt(roi, BLUR_KERNEL_LENGTH);
        
        % --- FIRST PASS: COARSE THRESHOLD --- %
        [roiBinarized, FIRST_PASS_HSV_MIN, FIRST_PASS_HSV_MAX, first_threshold] = firstPass(roiBlurred, ballColor, h, ...
            MAX_ITERATIONS, first_threshold, FIRST_PASS_HSV_MIN, FIRST_PASS_HSV_MAX, CLOSING_KERNEL_LENGTH, SMALLESTACCEPTABLEAREA);
        
       
        % --- CROP THE MULTIPLE ROI --- %
        [B, L] = bwboundaries(roiBinarized);
        stats = regionprops(L,'Centroid');
        
        if ~isempty(stats)
            for i=1:length(stats(:,1))
                % may be don't save candidateRoiRect 
                candidateRoiRect(i,:) = calcRoiSize(stats(i).Centroid, [size(frame,2), size(frame,1)], 4*BALL_SIZE, 2*BALL_SIZE);
                candidateBallsRoi{i} = imcrop(roiBlurred, candidateRoiRect(i, :));
            end
        else 
            candidateRoiRect = roiRect;
            candidateBallsRoi{1} = roiBlurred;
        end
        
        
        % --- SECOND PASS: RELAXED THRESHOLD --- % 
        [candidateBallsRoiBin,  SECOND_PASS_HSV_MIN, SECOND_PASS_HSV_MAX, second_threshold] = secondPass(candidateBallsRoi, ...
            ballColor, h, MAX_ITERATIONS, second_threshold, SECOND_PASS_HSV_MIN, SECOND_PASS_HSV_MAX, CLOSING_KERNEL_LENGTH, ...
            SMALLESTACCEPTABLEAREA, BALL_AREA);

        % --- FIND THE BALL! --- %
        for i = length(candidateBallsRoiBin)
            candidateRoiBin = candidateBallsRoiBin{i};

            % first try to find circles 
            % may be take this out, the ball is usually too small for this!
            [circles, radii] = getBallWithHough(candidateRoiBin, candidateRoiRect(i,:), RADIUSRANGE);

            if ~isempty(circles)
                if length(circles(:, 1)) == 1
                    ball_found = true;
                    if candidateRoiRect ~= roiRect
                        circles(1, 1) = circles(1,1) + roiRect(1);
                        circles(1, 2) = circles(1,2) + roiRect(2);
                    end
                    ball_position = circles; % radii?
                    %imshow(frame);
                    %viscircles(circles, radii);
                end
            end
            
            if ball_found == false 
                % second try
                posPositions = detectBallBoundaries(candidateRoiBin, candidateRoiRect(i,:));

                if ~isempty(posPositions)
                    if length(posPositions(:,1)) == 1
                        ball_found = true;
                        if candidateRoiRect ~= roiRect
                            posPositions(1, 1) = posPositions(1,1) + roiRect(1);
                            posPositions(1, 2) = posPositions(1,2) + roiRect(2);
                        end
                        ball_position = posPositions;
                    end
                end
            end

        end
        % RIGHT NOW THIS WILL GET THE LAST BALL FOUND; 
        % WILL BE CHANGED! 
        if ball_found == true
            positions = [positions ; ball_position];
        end
       
        frame_previous = frame;
        ball_found_prev = ball_found;        
    end
end

% let's plot the balls found, include the radius
hold on; plot(positions(:,1), positions(:,2), 'rx');

function roiRect = calcRoiSize(position, frame_size, ROI_WIDTH, ROI_HEIGHT)
    x_min = position(1) - (ROI_WIDTH / 2);
    if x_min < 0
        x_min = 0;
    end
    
    x_max = position(1) + (ROI_WIDTH / 2);
    if x_max > frame_size(1) - 1
        x_max = frame_size(1) - 1;
    end
    
    y_min = position(2) - (ROI_HEIGHT /2);
    if y_min < 0
        y_min = 0;
    end
    
    y_max = position(2) + (ROI_HEIGHT / 2);
    if y_max > frame_size(2) - 1
        y_max = frame_size(2) - 1;
    end
    
    roiRect = [x_min, y_min, x_max - x_min, y_max - y_min ];    
end

function roiBinarized = thresholdImg(roiBlurred, BALL_COLOR_HSV_MIN, BALL_COLOR_HSV_MAX, CLOSING_KERNEL_LENGTH, SMALLESTACCEPTABLEAREA)

    roi_hsv = rgb2hsv(roiBlurred);

    hueBand = roi_hsv (:,:,1);
    saturationBand = roi_hsv(:,:,2);
    intensityBand = roi_hsv(:,:,3); % intensity is in the whole range
    
    % Now apply each color band's particular thresholds to the color band
    hueMask = (hueBand >= BALL_COLOR_HSV_MIN(1)) & (hueBand <= BALL_COLOR_HSV_MAX(1));
    saturationMask = (saturationBand >= BALL_COLOR_HSV_MIN(2)) & (saturationBand <= BALL_COLOR_HSV_MAX(2));
    intensityMask = (intensityBand >= BALL_COLOR_HSV_MIN(3)) & (intensityBand <= BALL_COLOR_HSV_MAX(3));
    
    % Combine the masks to find where all 3 are "true."
    % Then we will have the mask of only the red parts of the image.
    ballColorMask = uint8(hueMask & intensityMask & saturationMask);
    
    % let's filter small objects
    %ballColorMask = uint8(bwareaopen(ballColorMask, SMALLESTACCEPTABLEAREA));
    
    % morphological closing operation
    structuringElement = strel('disk', CLOSING_KERNEL_LENGTH); 
    ballColorMask = imclose(ballColorMask, structuringElement);
    
    % Fill in any holes
    ballColorMask = uint8(imfill(ballColorMask, 'holes'));
    
    % convert the type of ballColorMask to the same data type as hueBand.
    roiBinarized = cast(ballColorMask, class(hueBand));
    
    % JUST FOR CODE COMPLETNESS, NOT USED IN THIS PROGRAM
    % Use the ball color object mask to mask out the ball color-only portions of the rgb image.
    %maskedImageH = ballColorMask .* hueBand;
    %maskedImageS = ballColorMask .* saturationBand;
    %maskedImageI = ballColorMask .* intensityBand;
    
    % Concatenate the masked color bands to form the rgb image.
    %roi_binarized = cat(3, maskedImageH, maskedImageS, maskedImageI);
    
    %imshow(roi_binarized);
end

function posPositions = detectBallBoundaries(roiBinarized, roiRect)
    [B, L] = bwboundaries(roiBinarized);
    %posPositions = [];
    %figure; imshow(frame);
    %hold on;
    posPositions = [];
    
    stats = regionprops(L,'Area','Centroid');

    threshold = 0.75; % let's see how to set this threshold better 

    % loop over the boundaries
    for k = 1:length(B)
        boundary = B{k};
        
        % compute a simple estimate of the object's perimeter
        delta_sq = diff(boundary).^2;
        perimeter = sum(sqrt(sum(delta_sq,2)));
        
        % obtain the area calculation corresponding to label 'k'
        area = stats(k).Area;
        
        % compute the roundness metric
        metric = 4*pi*area/perimeter^2;
        
        % display the results
        %metric_string = sprintf('%2.2f',metric);
        
        % mark objects above the threshold with a black circle
        if metric > threshold
            centroid = stats(k).Centroid;
            centroid(1) = centroid(1) + roiRect(1);
            centroid(2) = centroid(2) + roiRect(2);
            posPositions = [posPositions; centroid];
            %currentMax = metric;
            %plot(centroid(1),centroid(2),'ko');
        end

      %text(boundary(1,2)-35,boundary(1,1)+13,metric_string,'Color','y',...
       %    'FontSize',14,'FontWeight','bold');

    end

end

function [circles, radii] = getBallWithHough(roiBinarized, roiRect, RADIUSRANGE)
    cannyRoi = edge(roiBinarized, 'Canny');
    
    [circles, radii] = imfindcircles(cannyRoi, RADIUSRANGE); % dynamic set this threshold according to ball size!!
    
    % for each circel correct the center position 
   
    if ~isempty(circles)
        for i = 1:length(circles(:,1))
            circles(i, 1) = circles(i,1) + roiRect(1);
            circles(i, 2) = circles(i,2) + roiRect(2);
        end
    end
     
end

function [roiBinarized, BEST_FIRST_PASS_HSV_MIN, BEST_FIRST_PASS_HSV_MAX, first_threshold] = firstPass(roiBlurred, ballColor, h, ...
    MAX_ITERATIONS, first_threshold, FIRST_PASS_HSV_MIN, FIRST_PASS_HSV_MAX, CLOSING_KERNEL_LENGTH, SMALLESTACCEPTABLEAREA)

    iteration = 0;
    best_first_pass = [-1, 1000];
    candidate_found = 0;

    % FIRST PASS (COARSE THRESHOLD)
    % note for future: transform each pass in a function
    % question: should the treshold reset with each frame or use the
    % last one?
    while(~candidate_found && iteration < MAX_ITERATIONS)

        roiBinarized = thresholdImg(roiBlurred, FIRST_PASS_HSV_MIN, FIRST_PASS_HSV_MAX, CLOSING_KERNEL_LENGTH, SMALLESTACCEPTABLEAREA);

        candidateBallsInfo = bwconncomp(roiBinarized, 4);

        if candidateBallsInfo.NumObjects > 1
            if (best_first_pass(2) > candidateBallsInfo.NumObjects)
                best_first_pass = [iteration, candidateBallsInfo.NumObjects];
                BEST_FIRST_PASS_HSV_MIN =  FIRST_PASS_HSV_MIN;
                BEST_FIRST_PASS_HSV_MAX =  FIRST_PASS_HSV_MAX;
            end

            % increase the threshold
            first_threshold = round(first_threshold*0.8);

            if strcmp(ballColor, 'o')
                FIRST_PASS_HSV_MIN(1) = (h - first_threshold)/360;
                FIRST_PASS_HSV_MAX(1) = (h + first_threshold)/360;
            else
                %FIRST_PASS_HSV_MIN = [0, 0, 0.6];
                FIRST_PASS_HSV_MAX(2) = h+first_threshold/100;
            end

        elseif candidateBallsInfo.NumObjects == 1
            candidate_found = true;
            best_first_pass = [iteration,candidateBallsInfo.NumObjects];
            BEST_FIRST_PASS_HSV_MIN =  FIRST_PASS_HSV_MIN;
            BEST_FIRST_PASS_HSV_MAX =  FIRST_PASS_HSV_MAX;
        else
            % decrease the threshold
            first_threshold = round(first_threshold * 1.2);
            if strcmp(ballColor, 'o')
                FIRST_PASS_HSV_MIN(1) = (h - first_threshold)/360;
                FIRST_PASS_HSV_MAX(1) = (h + first_threshold)/360;

                % handle edges
                if(FIRST_PASS_HSV_MIN(1) < 0)
                    FIRST_PASS_HSV_MIN(1) = 0;
                end
                if(FIRST_PASS_HSV_MAX(1) > 1)
                    FIRST_PASS_HSV_MAX(1) = 1;
                end
            else
                %FIRST_PASS_HSV_MIN = [0, 0, 0.60];
                FIRST_PASS_HSV_MAX(2) = h+first_threshold/100;
                if(FIRST_PASS_HSV_MAX(1) > 1)
                    FIRST_PASS_HSV_MAX(1) = 1;
                end
            end
        end
        iteration = iteration + 1;
    end

    % RESTORE BEST_THRESHOLD

    if (~candidate_found && iteration - 1 ~= best_first_pass(1))
        if best_first_pass(2) < 1000
            roiBinarized = thresholdImg(roiBlurred, BEST_FIRST_PASS_HSV_MIN, BEST_FIRST_PASS_HSV_MAX, CLOSING_KERNEL_LENGTH, SMALLESTACCEPTABLEAREA);
        else 
            BEST_FIRST_PASS_HSV_MIN = FIRST_PASS_HSV_MIN;
            BEST_FIRST_PASS_HSV_MAX = FIRST_PASS_HSV_MAX;
        end
    end

end

function [candidateBallsRoi,  BEST_SECOND_PASS_HSV_MIN, BEST_SECOND_PASS_HSV_MAX, second_threshold] = secondPass(candidateBallsRoi, ...
   ballColor, h, MAX_ITERATIONS, second_threshold, SECOND_PASS_HSV_MIN, SECOND_PASS_HSV_MAX, CLOSING_KERNEL_LENGTH, ...
   SMALLESTACCEPTABLEAREA, BALL_AREA )

    % SECOND PASS AND CHOOSE BALL
    for i = length(candidateBallsRoi)


        % SECOND PASS (RELAXED THRESHOLD)
        iteration = 0;
        area_ok = false;
        candidateRoi = candidateBallsRoi{i};
        best_sec_pass = [-1, Inf];

        while(~area_ok && iteration < MAX_ITERATIONS)

            candidateRoiBinarized = thresholdImg(candidateRoi, SECOND_PASS_HSV_MIN, SECOND_PASS_HSV_MAX, CLOSING_KERNEL_LENGTH, SMALLESTACCEPTABLEAREA);

            [~, L] = bwboundaries(candidateRoiBinarized, 'noholes');
            stats = regionprops(L,'Area');
            
                        
            if ~isempty(stats)
                % confirm if this is okay
                area = max([stats.Area]);
                aux_area = (BALL_AREA - area)/BALL_AREA;
                % confirm this area threshold
                if  (aux_area < 0.4) && (aux_area > -0.4)
                    area_ok = true;
                    BEST_SECOND_PASS_HSV_MAX = SECOND_PASS_HSV_MAX;
                    BEST_SECOND_PASS_HSV_MIN = SECOND_PASS_HSV_MIN;

                    % untested control sequence
                elseif area > BALL_AREA

                    if best_sec_pass(2) > area - BALL_AREA
                        best_sec_pass = [iteration, area - BALL_AREA];
                        BEST_SECOND_PASS_HSV_MAX = SECOND_PASS_HSV_MAX;
                        BEST_SECOND_PASS_HSV_MIN = SECOND_PASS_HSV_MIN;
                    end

                    % increase the threshold
                    second_threshold = round(second_threshold*0.85);

                    if strcmp(ballColor, 'o')
                        SECOND_PASS_HSV_MIN(1) = (h - second_threshold)/360;
                        SECOND_PASS_HSV_MAX(1) = (h + second_threshold)/360;
                    else
                        SECOND_PASS_HSV_MAX(2) = h+second_threshold/100;
                    end
                else

                    if best_sec_pass(2) > BALL_AREA - area
                        best_sec_pass = [iteration, BALL_AREA - area];
                        BEST_SECOND_PASS_HSV_MAX = SECOND_PASS_HSV_MAX;
                        BEST_SECOND_PASS_HSV_MIN = SECOND_PASS_HSV_MIN;
                    end

                    % decrease the threshold
                    second_threshold = round(second_threshold * 1.15);
                    if strcmp(ballColor, 'o')
                        SECOND_PASS_HSV_MIN(1) = (h - second_threshold)/360;
                        SECOND_PASS_HSV_MAX(1) = (h + second_threshold)/360;

                        % handle edges
                        if(SECOND_PASS_HSV_MIN(1) < 0)
                            SECOND_PASS_HSV_MIN(1) = 0;
                        end
                        if(SECOND_PASS_HSV_MAX(1) > 1)
                            SECOND_PASS_HSV_MAX(1) = 1;
                        end
                    else
                        SECOND_PASS_HSV_MAX(2) = h+second_threshold/100;
                        if(SECOND_PASS_HSV_MAX(1) > 1)
                            SECOND_PASS_HSV_MAX(1) = 1;
                        end
                    end
                end
            else
                % nothing found, decrease the threshold
                second_threshold = round(second_threshold * 1.15);
                if strcmp(ballColor, 'o')
                    SECOND_PASS_HSV_MIN(1) = (h - second_threshold)/360;
                    SECOND_PASS_HSV_MAX(1) = (h + second_threshold)/360;
                    
                    % handle edges
                    if(SECOND_PASS_HSV_MIN(1) < 0)
                        SECOND_PASS_HSV_MIN(1) = 0;
                    end
                    if(SECOND_PASS_HSV_MAX(1) > 1)
                        SECOND_PASS_HSV_MAX(1) = 1;
                    end
                else
                    SECOND_PASS_HSV_MAX(2) = h+second_threshold/100;
                    if(SECOND_PASS_HSV_MAX(1) > 1)
                        SECOND_PASS_HSV_MAX(1) = 1;
                    end
                end
                
                
            end
                
            iteration = iteration + 1;
        end

        % RESTORE BEST THRESHOLD
        if (~area_ok && iteration - 1 ~= best_sec_pass(1))
            if best_sec_pass(2) < Inf
                candidateRoiBinarized = thresholdImg(candidateRoi, BEST_SECOND_PASS_HSV_MIN, BEST_SECOND_PASS_HSV_MAX, CLOSING_KERNEL_LENGTH, SMALLESTACCEPTABLEAREA);
            else
                BEST_SECOND_PASS_HSV_MIN = SECOND_PASS_HSV_MIN;
                BEST_SECOND_PASS_HSV_MAX = SECOND_PASS_HSV_MAX;
            end
        end

        candidateBallsRoi{i} = candidateRoiBinarized;
    end


end
