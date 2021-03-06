clc;
clear;
close all;
fontSize = 20;


MAX_ITERATIONS = 3;

SMALLESTACCEPTABLEAREA = 40;
BALL_SIZE = 20; % pixels of the ball, we should determine this...
RADIUSRANGE = [round(BALL_SIZE/2-0.3*BALL_SIZE/2), round(BALL_SIZE/2+0.3*BALL_SIZE/2)]; % PAY ATTENTION, it's THE RADIUS!
ROI_WIDTH = BALL_SIZE * 10; 
ROI_HEIGHT = BALL_SIZE * 5;

% kernelfor blurring. VALUE: 3-sigma rule
BLUR_KERNEL_LENGTH = (BALL_SIZE-1)/6; 

CLOSING_KERNEL_LENGTH = round(BALL_SIZE/2); % study a little more this..

% adapt this to the color of the ball
% orange
h = 30; h_threshold = 20;
BALL_COLOR_HSV_MIN = [(h - h_threshold)/180, 50/255, 50/255];
BALL_COLOR_HSV_MAX = [(h + h_threshold)/180, 255/255, 255/255];

% white
%h = 30; h_threshold = 20;
%BALL_COLOR_HSV_MIN = [0, 0/255, 0.70];
%BALL_COLOR_HSV_MAX = [1, h+h_treshold/255, 1];

%prompt = 'Input the name of the video file: ';
%videoName = input(prompt, 's');

videoName = 'Video/orange_ball_high.mp4';

v = VideoReader(videoName);

%prompt = 'Frames to be skipped: ';
%i = input(prompt, 's');

skips = 5;
i = 0;

positions = []; 
ball_found_prev = false; 

first_frame = true;

while hasFrame(v)
    
    % skip the bad frames at the beginning
    if i < skips  
        if hasFrame(v)
            frame = readFrame(v);
        end
        i=i+1;
    else
        
        frame = readFrame(v);
        ball_found = 0;
        roiRect = [0, 0, 0, 0];
        
        % calculate the ROI 
        if ball_found_prev == 1
            roiRect = calcRoiSize(positions(length(positions(:,1)), :), [size(frame,2), size(frame,1)], ROI_WIDTH, ROI_HEIGHT);
            roi = imcrop(frame, roiRect);
        else
            % we ask the user for the initial location of the ball
            % helps in low res cases
            if first_frame  == true
                imshow(frame);
                [pos(1), pos(2)] = getpts;
                roiRect = calcRoiSize(pos, [size(frame,2), size(frame,1)], ROI_WIDTH, ROI_HEIGHT);
                roi = imcrop(frame, roiRect);
                first_frame = false;
                
            % this needs to be enhenced 
            else         
               roi = frame;
            end
        end
          
        iteration = 0;
        roiBlurred = imgaussfilt(roi, BLUR_KERNEL_LENGTH);
        best_iteration = [0, 1000];
        
        % first pass 
        while(~ball_found && iteration < MAX_ITERATIONS)
        
            roiBinarized = thresholdImg(roiBlurred, BALL_COLOR_HSV_MIN, BALL_COLOR_HSV_MAX, CLOSING_KERNEL_LENGTH, SMALLESTACCEPTABLEAREA);

            % first try to find circles 
            [circles, radii] = getCandidateBalls(roiBinarized, roiRect, RADIUSRANGE);

            if ~isempty(circles)
                if length(circles(:, 1)) == 1
                    ball_found = true;
                    ball_position = circles;
                    %imshow(frame);
                    %viscircles(circles, radii);
                end
                
                if (best_iteration(2) > length(circles(:,1)))
                    best_iteration = [i, length(circles(:,1))];
                    BEST_BALL_COLOR_HSV_MIN =  BALL_COLOR_HSV_MIN;
                    BEST_BALL_COLOR_HSV_MAX =  BALL_COLOR_HSV_MAX;
                end
            end
            
            if ball_found == false 
                % second try
                posPositions = detectBallBoundaries(roiBinarized, roiRect);

                if ~isempty(posPositions)
                    if length(posPositions(:,1)) == 1
                        ball_found = true;
                        ball_position = posPositions;
                    else 
                        h_threshold = round(h_threshold*0.8);
                        BALL_COLOR_HSV_MIN = [(h - h_threshold)/180, 50/255, 50/255];
                        BALL_COLOR_HSV_MAX = [(h + h_threshold)/180, 255/255, 255/255];
                    end
                    
                    if (best_iteration(2) > length(posPositions(:,1)) && ~isempty(posPositions(:,1)))
                        best_iteration = [i, length(posPositions(:,1))];
                        BEST_BALL_COLOR_HSV_MIN =  BALL_COLOR_HSV_MIN;
                        BEST_BALL_COLOR_HSV_MAX =  BALL_COLOR_HSV_MAX;
                    end
                else
                    h_threshold = round(h_threshold * 1.2);
                    BALL_COLOR_HSV_MIN = [(h - h_threshold)/180, 50/255, 50/255];
                    if(BALL_COLOR_HSV_MIN(1) < 0)
                            BALL_COLOR_HSV_MIN(1) = 0;
                    end
                    BALL_COLOR_HSV_MAX = [(h + h_threshold)/180, 255/255, 255/255];
                    if(BALL_COLOR_HSV_MAX(1) > 1)
                        BALL_COLOR_HSV_MAX(1) = 1;
                    end
                end
            end
            
            iteration = iteration + 1;
        end 
        
        
%         if (~ball_found && iteration + 1 ~= best_iteration(1))
%             BALL_COLOR_HSV_MIN =  BEST_BALL_COLOR_HSV_MIN;
%             BALL_COLOR_HSV_MAX =  BEST_BALL_COLOR_HSV_MAX;
%             
%         end
        
        
        if ball_found == true
            positions = [positions ; ball_position];
        end
        
        % here we start displaying the frames.. 
        % not implemented yet
        %if(ball_found_prev)
        %    drawTrackingWithPrevious(frame, positions); 
        %else
        %    drawTrackingWithoutPrevious(frame, positions, roiRect);
        %end
        
        frame_previous = frame;
        ball_found_prev = ball_found;
        
    end
end

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

function roiBinarized = thresholdImg(inputRoi, BALL_COLOR_HSV_MIN, BALL_COLOR_HSV_MAX, CLOSING_KERNEL_LENGTH, SMALLESTACCEPTABLEAREA)

    roi_hsv = rgb2hsv(inputRoi);

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
    ballColorMask = uint8(bwareaopen(ballColorMask, SMALLESTACCEPTABLEAREA));
    
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

function [circles, radii] = getCandidateBalls(roiBinarized, roiRect, RADIUSRANGE)
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

function drawTrackingWithPrevious(frame, positions)
        
    figure;
    imshow(frame);
    hold on;
    
    if size(positions, 1) > 1
        for i = 1:(length(positions(:,1))-1)    
            plot( [positions(i, 1), positions(i+1, 1)], [positions(i, 2), positions(i+1, 2)]);
        end
    end
    
end

function drawTrackingWithoutPrevious(frame, positions, roiRect)
        
    drawTrackingWithPrevious(frame, positions);
    rectangle('Position', roiRect);
    
end

