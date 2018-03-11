clc;
clear;
close all;
fontSize = 20;

MAX_ITERATIONS = 3;

% ---- BALL RELATED CONSTANTS ---- % 

BALL_SIZE = 13; % pixels of the ball, CALIBRATE THE CAMERA, this comes automatically 
BALL_AREA = pi * (BALL_SIZE/2)^2;
BALL_PERIMETER = 2 * pi * (BALL_SIZE/2);
RADIUSRANGE = [round(BALL_SIZE/2-0.3*BALL_SIZE/2), round(BALL_SIZE/2+0.3*BALL_SIZE/2)]; 
ROI_WIDTH = BALL_SIZE * 10; %this can vary with the video FPS 
ROI_HEIGHT = BALL_SIZE * 5;
ADAPTIVE_WIDTH = ROI_WIDTH;
ADAPTIVE_HEIGHT = ROI_HEIGHT;

% kernel for blurring. VALUE: 3-sigma rule
BLUR_KERNEL_LENGTH = (BALL_SIZE-1)/6; 
CLOSING_KERNEL_LENGTH = round(BALL_SIZE/2);

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
    h = 35; first_threshold = 10; second_threshold = 25;
    FIRST_PASS_HSV_MIN = [(h - first_threshold)/360, 45/100, 30/100];
    FIRST_PASS_HSV_MAX = [(h + first_threshold)/360, 1, 1];
    % may be add a variable in the saturation as well 
    SECOND_PASS_HSV_MIN = [(h - second_threshold)/360, 30/100, 0];
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
                [positions(1), positions(2)] = getpts;
                positions(3) = v.CurrentTime;
                roiRect = calcRoiSize(positions, [size(frame,2), size(frame,1)], ROI_WIDTH, ROI_HEIGHT);
                frame = readFrame(v);
                roi = imcrop(frame, roiRect);
                first_frame = false;
                
            else
               % MAY BE IMPLEMENT linear EXTRAPOLATION!!
               ADAPTIVE_WIDTH = 1.8 * ADAPTIVE_WIDTH;
               ADAPTIVE_HEIGHT = 1.8 * ADAPTIVE_HEIGHT;
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
            MAX_ITERATIONS, first_threshold, FIRST_PASS_HSV_MIN, FIRST_PASS_HSV_MAX, CLOSING_KERNEL_LENGTH);
       
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
        [candidateBallsRoiBin, candidateBallsHSV] = secondPass(candidateBallsRoi, ...
            ballColor, h, MAX_ITERATIONS, second_threshold, SECOND_PASS_HSV_MIN, SECOND_PASS_HSV_MAX, CLOSING_KERNEL_LENGTH, ...
            BALL_AREA);
        
        posError = [0, 0, inf];
        
        % --- FIND THE BALL! --- %
        for i = 1:length(candidateBallsRoiBin)
            candidateRoiBin = candidateBallsRoiBin{i};
            
            % in case we're using a ROI inside the ROI
            if candidateRoiRect ~= roiRect
                candidateRoiRect(i, 1) = candidateRoiRect(i, 1) + roiRect(1);
                candidateRoiRect(i, 2) = candidateRoiRect(i, 2) + roiRect(2);
            end
            
            nextTs = v.CurrentTime; %+1/v.FrameRate;
            error = twoStageBallDetection(candidateRoiBin, candidateRoiRect(i,:), positions, nextTs, BALL_SIZE, BALL_AREA, BALL_PERIMETER);
            
            
            if error < posError(3)
                posError = error;
                SECOND_PASS_HSV_MIN = candidateBallsHSV(i, 1:3);
                SECOND_PASS_HSV_MAX = candidateBallsHSV(i, 4:6);
                second_threshold = candidateBallsHSV(i, 7);
            end

        end
        
        if posError(3) < inf 
            ball_found = true;
            positions = [positions ; posError(1), posError(2), v.CurrentTime];
        end
         
        frame_previous = frame;
        ball_found_prev = ball_found;        
    end
end

% let's plot the balls found, include the radius
%hold on; plot(positions(:,1), positions(:,2), 'rx');
interval = positions(1,3):0.0002:positions(length(positions),3);
vq = interp1(positions(:, 3), positions(:, 1:2), interval , 'spline');
hold on; plot(vq(:,1), vq(:,2), 'r');

% search for the bounce 
% analyzing y in 5 frames
candidateTs = [];
for i= 3:length(positions)-2
   if ((positions(i-2, 2) < positions(i-1, 2)) && (positions(i-1, 2) < positions(i,2)) && (positions(i, 2) > positions(i+1, 2)) && (positions(i+1, 2) > positions(i+2,2))) 
       candidateTs = [candidateTs; i];
   end
end

% estimate the inside the frames point of bounce 
bounceTs = [];
bouncePoint = [];
for j = 1:length(candidateTs)
    interval = positions(candidateTs(j)-1,3):0.0002:positions(candidateTs(j)+1,3);
    vq = interp1(positions(:, 3), positions(:, 1:2), interval , 'spline');
    for i = 4:length(vq)-3
       if ((vq(i-3, 2) < vq(i-2, 2)) && (vq(i-2, 2) < vq(i-1, 2)) && (vq(i-1, 2) < vq(i,2)) && (vq(i, 2) > vq(i+1, 2)) && (vq(i+1, 2) > vq(i+2,2)) && (vq(i+2, 2) > vq(i+3,2))) 
           bounceTs = [bounceTs ; interval(i)];
           bouncePoint = [bouncePoint; vq(i, :)];
       end
    end
end


% ---- END OF MAIN LOOP ---- % 

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

% still have to look the SMALLESTACCEPTABLEAREA issue 
function roiBinarized = thresholdImg(roiBlurred, BALL_COLOR_HSV_MIN, BALL_COLOR_HSV_MAX, CLOSING_KERNEL_LENGTH)

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

function [xc,yc,R,a] = circlefit(x,y)
%CIRCFIT  Fits a circle in x,y plane
%
% [XC, YC, R, A] = CIRCFIT(X,Y)
% Result is center point (yc,xc) and radius R.  A is an optional
% output describing the circle's equation:
%
%   x^2+y^2+a(1)*x+a(2)*y+a(3)=0
         
% by Bucher izhak 25/oct/1991

n=length(x);  xx=x.*x; yy=y.*y; xy=x.*y;
A=[sum(x) sum(y) n;sum(xy) sum(yy) sum(y);sum(xx) sum(xy) sum(x)];
B=[-sum(xx+yy) ; -sum(xx.*y+yy.*y) ; -sum(xx.*x+xy.*y)];
a=A\B;
xc = -.5*a(1);
yc = -.5*a(2);
R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));

end

function minError = twoStageBallDetection(roiBinarized, roiRect, positions, nextTs, BALL_SIZE, BALL_AREA, BALL_PERIMETER)
    % weights for the second stage equation 
    ww = 0.2;
    wa = 0.2;
    wh = 0.2;
    wp = 0.2;
    wr = 0.2;

    % bwboundaries gives ROWS 3and COLUMNS = Height, width
    [B, L] = bwboundaries(roiBinarized);
    minError = [0, 0, inf]; 
    stats = regionprops(L, 'Area', 'Centroid', 'BoundingBox');
   
    % loop over the boundaries
    for k = 1:length(B)
        
        boundary = B{k};
        
        % ----- FIRST STAGE ----- %
        
        % --- ROUNDED UPPER CONTOUR --- %
        % Get upper contour pixels ( y = boundary(:,1) )
        Eruc = 0;
        upperContour = boundary( boundary(:,1) <= stats(k).BoundingBox(4)/2 + stats(k).BoundingBox(2), :);
        % fit in a circle
        if length(upperContour) < 4
            RUC = false;
        else
            [xc, yc, R] = circlefit(upperContour(:,2), upperContour(:,1));
            
            % calculate the error
            n = length(upperContour(:,1));
            sumError = 0;
            for i=1:n
                d = sqrt((upperContour(i,2)- xc)^2 + (upperContour(i,1) - yc)^2);
                sumError = sumError + abs((d - R)/R);
            end
            
            Eruc = sumError/n;
            % fine tune the threshold
            if Eruc < 0.1
                RUC = true;
            else
                RUC = false;
            end
        end
        
        % pay attention how this will work in practice, there may be
        % problem until a prediciont can be made!
        centroid = stats(k).Centroid;
        centroid(1) = centroid(1) + roiRect(1);
        centroid(2) = centroid(2) + roiRect(2);
        
        mc = false;
        mp = false;
        T = false;
        if ~isempty(positions)
            if length(positions(:,1)) > 1
                
                prediction = interp1(positions(:,3), positions(:,1:2), nextTs, 'linear', 'extrap');
                
                % --- POSITION --- %
                distOOI = pdist([centroid; prediction], 'euclidean');
                if distOOI < 1.5*BALL_SIZE % fine tune the maximum dist
                    T = true;
                end
                
                % --- MOTION --- %
                % taking into account the last known ball location..
                distC = pdist([centroid; positions(length(positions(:,1)), 1:2)], 'euclidean');
                if distC > 0.5 * BALL_SIZE % set this threshold!
                    mc = true;
                end
                
                %distP = pdist([prediction; positions(length(positions(:,1)), 1:2)], 'euclidean');
                %if distP > 0.5 * BALL_SIZE
                %    mp = true;
                %end
            elseif length(positions(:,1)) == 1
                % still don't have enough for a prediction, check only
                % movement
                distC = pdist([centroid; positions(1,1:2)], 'euclidean');
                if distC > 0.5 * BALL_SIZE % set this threshold!
                    mc = true;
                end
                %let's give a free point
                mp = 1;
            end
        end
        
        % ----- SECOND STAGE ----- %
        if RUC + T + mc + mp > 1
            
            area = stats(k).Area;
            
            % compute a simple estimate of the object's perimeter
            delta_sq = diff(boundary).^2;
            perimeter = sum(sqrt(sum(delta_sq,2)));
            
            % roundness metric
            roundness = 4*pi*area/perimeter^2;
            
            
            maxHeight = stats(k).BoundingBox(4);
            maxWidth = stats(k).BoundingBox(3);
            
            E = (Eruc + wa * abs(area - BALL_AREA)/BALL_AREA + ww * abs(maxWidth - BALL_SIZE)/(BALL_SIZE) + ...
                wh * abs(maxHeight - BALL_SIZE)/(BALL_SIZE) + wp * abs(perimeter - BALL_PERIMETER)/BALL_PERIMETER + ...
                wr * abs(roundness - 1))/ (5+1)*(RUC + T + mc + mp);
            
            if (E < minError(3) && E < 0.48)
                minError(3) = E;
                minError(1) = centroid(1);
                minError(2) = centroid(2);
            end
        end
        
    end
    
    % should return the error function of the object 
end

function [roiBinarized, BEST_FIRST_PASS_HSV_MIN, BEST_FIRST_PASS_HSV_MAX, best_threshold] = firstPass(roiBlurred, ballColor, h, ...
    MAX_ITERATIONS, first_threshold, FIRST_PASS_HSV_MIN, FIRST_PASS_HSV_MAX, CLOSING_KERNEL_LENGTH)

    iteration = 0;
    best_first_pass = [-1, 1000];
    candidate_found = 0;

    % FIRST PASS (COARSE THRESHOLD)
    % note for future: transform each pass in a function
    % question: should the treshold reset with each frame or use the
    % last one?
    while(~candidate_found && iteration < MAX_ITERATIONS)

        roiBinarized = thresholdImg(roiBlurred, FIRST_PASS_HSV_MIN, FIRST_PASS_HSV_MAX, CLOSING_KERNEL_LENGTH);

        candidateBallsInfo = bwconncomp(roiBinarized, 4);

        if candidateBallsInfo.NumObjects > 1
            if (best_first_pass(2) >= candidateBallsInfo.NumObjects)
                best_first_pass = [iteration, candidateBallsInfo.NumObjects];
                BEST_FIRST_PASS_HSV_MIN =  FIRST_PASS_HSV_MIN;
                BEST_FIRST_PASS_HSV_MAX =  FIRST_PASS_HSV_MAX;
                best_threshold = first_threshold;
            end

            % increase the threshold
            first_threshold = first_threshold*0.80;

            if strcmp(ballColor, 'o')
                FIRST_PASS_HSV_MIN(1) = (h - first_threshold)/360;
                %FIRST_PASS_HSV_MIN(2) = (75 - first_threshold/2)/100;
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
            best_threshold = first_threshold;
        else
            % decrease the threshold
            first_threshold = first_threshold * 1.2;
            if strcmp(ballColor, 'o')
                FIRST_PASS_HSV_MIN(1) = (h - first_threshold)/360;
                %FIRST_PASS_HSV_MIN(2) = (75 - first_threshold/2)/100;
                FIRST_PASS_HSV_MAX(1) = (h + first_threshold)/360;

                % handle edges
                if(FIRST_PASS_HSV_MIN(1) < 0)
                    FIRST_PASS_HSV_MIN(1) = 0;
                end
                %if FIRST_PASS_HSV_MIN(2) < 0
                %    FIRST_PASS_HSV_MIN(2) = 0;
                %end
                if(FIRST_PASS_HSV_MAX(1) > 1)
                    FIRST_PASS_HSV_MAX(1) = 1;
                end
            else
                %FIRST_PASS_HSV_MIN = [0, 0, 0.60];
                FIRST_PASS_HSV_MAX(2) = h+first_threshold/100;
                if(FIRST_PASS_HSV_MAX(2) > 1)
                    FIRST_PASS_HSV_MAX(2) = 1;
                end
            end
        end
        iteration = iteration + 1;
    end

    % RESTORE BEST_THRESHOLD

    if (~candidate_found && iteration - 1 ~= best_first_pass(1))
        if best_first_pass(2) < 1000
            roiBinarized = thresholdImg(roiBlurred, BEST_FIRST_PASS_HSV_MIN, BEST_FIRST_PASS_HSV_MAX, CLOSING_KERNEL_LENGTH);
        else 
            % we could give additional tries here with a smaller
            % threshold.. 
            BEST_FIRST_PASS_HSV_MIN = FIRST_PASS_HSV_MIN;
            BEST_FIRST_PASS_HSV_MAX = FIRST_PASS_HSV_MAX;
            best_threshold = first_threshold;
        end
    end

end

function [candidateBallsRoi, candidateBallsHSV] = secondPass(candidateBallsRoi, ...
   ballColor, h, MAX_ITERATIONS, second_threshold, SECOND_PASS_HSV_MIN, SECOND_PASS_HSV_MAX, CLOSING_KERNEL_LENGTH, ...
   BALL_AREA )

    %best_sec_pass = Inf; 
    candidateBallsHSV = [];

    % SECOND PASS (RELAXED THRESHOLD)  
    for i = 1:length(candidateBallsRoi)

        iteration = 0;
        area_ok = false;
        candidateRoi = candidateBallsRoi{i};
        %best sec pass = interation, threshold
        best_it = [-1, Inf];
        IT_HSV_MIN = SECOND_PASS_HSV_MIN;
        IT_HSV_MAX = SECOND_PASS_HSV_MAX;
        it_threshold = second_threshold;
        
        while(~area_ok && iteration < MAX_ITERATIONS)

            candidateRoiBinarized = thresholdImg(candidateRoi, IT_HSV_MIN, IT_HSV_MAX, CLOSING_KERNEL_LENGTH);

            [~, L] = bwboundaries(candidateRoiBinarized, 'noholes');
            stats = regionprops(L,'Area');
            
            if ~isempty(stats)
                % not ok!
                area = max([stats.Area]);
                aux_area = abs(BALL_AREA - area)/BALL_AREA;
                % confirm this area threshold
                if  aux_area < 0.4
                    area_ok = true;
                    BEST_ITERATION_HSV_MAX = IT_HSV_MAX;
                    BEST_ITERATION_HSV_MIN = IT_HSV_MIN;
                    best_iteration_threshold = it_threshold;
                    best_it = [iteration, abs(area - BALL_AREA)];
                    % untested control sequence
                elseif area > BALL_AREA

                    if best_it(2) > area - BALL_AREA
                        best_it = [iteration, area - BALL_AREA];
                        BEST_ITERATION_HSV_MAX = IT_HSV_MAX;
                        BEST_ITERATION_HSV_MIN = IT_HSV_MIN;
                        best_iteration_threshold = it_threshold;
                    end

                    % increase the threshold
                    it_threshold = it_threshold*0.80;

                    if strcmp(ballColor, 'o')
                        IT_HSV_MIN(1) = (h - it_threshold)/360;
                        %SECOND_PASS_HSV_MIN(2) = (55 - second_threshold/2)/100;
                        IT_HSV_MAX(1) = (h + it_threshold)/360;
                    else
                        IT_HSV_MAX(2) = h + it_threshold/100;
                    end
                else

                    if best_it(2) > BALL_AREA - area
                        best_it = [iteration, BALL_AREA - area];
                        BEST_ITERATION_HSV_MAX = IT_HSV_MAX;
                        BEST_ITERATION_HSV_MIN = IT_HSV_MIN;
                        best_iteration_threshold = it_threshold;
                    end

                    % decrease the threshold
                    it_threshold = it_threshold * 1.2;
                    if strcmp(ballColor, 'o')
                        IT_HSV_MIN(1) = (h - it_threshold)/360;
                        %SECOND_PASS_HSV_MIN(2) = (55 - second_threshold/2)/100;
                        IT_HSV_MAX(1) = (h + it_threshold)/360;

                        % handle edges
                        if(IT_HSV_MIN(1) < 0)
                            IT_HSV_MIN(1) = 0;
                        end
                        %if SECOND_PASS_HSV_MIN(2) < 0
                        %    SECOND_PASS_HSV_MIN(2) = 0;
                        %end
                        if(IT_HSV_MAX(1) > 1)
                            IT_HSV_MAX(1) = 1;
                        end
                    else
                        IT_HSV_MAX(2) = h+it_threshold/100;
                        if(IT_HSV_MAX(1) > 1)
                            IT_HSV_MAX(1) = 1;
                        end
                    end
                end
            else
                % nothing found, decrease the threshold
                it_threshold = it_threshold * 1.2;
                if strcmp(ballColor, 'o')
                    IT_HSV_MIN(1) = (h - it_threshold)/360;
                    %SECOND_PASS_HSV_MIN(2) = (55 - second_threshold/2)/100;
                    IT_HSV_MAX(1) = (h + it_threshold)/360;
                    
                    % handle edges
                    if(IT_HSV_MIN(1) < 0)
                        IT_HSV_MIN(1) = 0;
                    end
                    %if SECOND_PASS_HSV_MIN(2) < 0
                    %    SECOND_PASS_HSV_MIN(2) = 0;
                    %end
                    if(IT_HSV_MAX(1) > 1)
                        IT_HSV_MAX(1) = 1;
                    end
                else
                    IT_HSV_MAX(2) = h+it_threshold/100;
                    if(IT_HSV_MAX(1) > 1)
                        IT_HSV_MAX(1) = 1;
                    end
                end
                
                
            end
                
            iteration = iteration + 1;
        end

        % RESTORE BEST THRESHOLD
        if (~area_ok && iteration - 1 ~= best_it(1))
            if best_it(2) < Inf
                candidateRoiBinarized = thresholdImg(candidateRoi, BEST_ITERATION_HSV_MIN, BEST_ITERATION_HSV_MAX, CLOSING_KERNEL_LENGTH);
            else
                BEST_ITERATION_HSV_MIN = IT_HSV_MIN;
                BEST_ITERATION_HSV_MAX = IT_HSV_MAX;
                best_iteration_threshold = it_threshold;
            end
        end
       
%         if best_sec_pass > best_it(2)
%             best_sec_pass = best_it(2);
%             BEST_SEC_PASS_HSV_MAX = BEST_ITERATION_HSV_MAX;
%             BEST_SEC_PASS_HSV_MIN = BEST_ITERATION_HSV_MIN;
%             best_sec_threshold = best_iteration_threshold;
%         elseif best_it(2) == Inf
%             BEST_SEC_PASS_HSV_MAX = SECOND_PASS_HSV_MAX;
%             BEST_SEC_PASS_HSV_MIN = SECOND_PASS_HSV_MIN;
%             best_sec_threshold = second_threshold;
%         end
         candidateBallsRoi{i} = candidateRoiBinarized;
        % EITHER KEEP BEST_IT OR THIS BELOW!! 
        if best_it(2) == inf
            candidateBallsHSV = [candidateBallsHSV; SECOND_PASS_HSV_MIN, SECOND_PASS_HSV_MAX, second_threshold];
        else
            candidateBallsHSV = [candidateBallsHSV; BEST_ITERATION_HSV_MIN, BEST_ITERATION_HSV_MAX, best_iteration_threshold];
        end
    end

end
