clc;
clear;
close all;
fontSize = 20;

SMALLESTACCEPTABLEAREA = 20;
BALL_SIZE = 20; % pixels of the ball, we should determine this...
RADIUSRANGE = [round(BALL_SIZE/2-0.3*BALL_SIZE/2), round(BALL_SIZE/2+0.3*BALL_SIZE/2)]; % PAY ATTENTION, it's THE RADIUS!
ROI_WIDTH = BALL_SIZE * 15; 
ROI_HEIGHT = BALL_SIZE * 7;


% kernels for blurring 
BLUR_KERNEL_LENGTH = 2;
CLOSING_KERNEL_LENGTH = 27;

% adapt this to the color of the ball
% if orange (RGB = [255, 252, 31]), HSV is = [30, 224, 255]
% if white (RGB = [255, 255, 255]), HSV is = [0, 0, 255]
% ignore V of the HSV to not depend on the light intensity
% orange
h = 30; h_threshold = 30;
BALL_COLOR_HSV_MIN = [(h - h_threshold)/180, 50/255, 0];
BALL_COLOR_HSV_MAX = [(h + h_threshold)/180, 255/255, 1];

% white
%h = 0; h_threshold = 20;
%BALL_COLOR_HSV_MIN = [0, 0/255, 0.70];
%BALL_COLOR_HSV_MAX = [1, 60/255, 1];


%prompt = 'Input the name of the video file: ';
%videoName = input(prompt, 's');

frame = imread('orange_ball_clip.jpg');

figure;
imshow(frame);
hold on;

%positions = []; 
%ball_found_prev = false; 

[positions(1), positions(2)] = getpts();
ball_found_prev = true;

ball_found = 0;
roiRect = [0, 0, 0, 0];

% we may need to resize frames

% if the ball was previously found, resize ROI, else the whole
% frame
if ball_found_prev == 1
    roiRect = calcRoiSize(positions(length(positions(:,1)), :), [size(frame, 2), size(frame, 1)], ROI_WIDTH, ROI_HEIGHT);
    roi = imcrop(frame, roiRect);
else
    roi = frame;
end

% let's study if this blurring is necessary
roiBlurred = imgaussfilt(roi, BLUR_KERNEL_LENGTH);

% first try
roiBinarized = thresholdImg(roiBlurred, BALL_COLOR_HSV_MIN, BALL_COLOR_HSV_MAX, CLOSING_KERNEL_LENGTH, SMALLESTACCEPTABLEAREA);
[circles, radii] = getCandidateBalls(roiBinarized, roiRect, RADIUSRANGE);

% if only one circle, that's the one
if ~isempty(circles)
    if length(circles(:, 1)) == 1
        ball_found = true;
        ball_position = circles;
        imshow(frame);
        viscircles(circles, radii);
    end
end

if ball_found == false
    % second try
    %the ball is not round enough (or other balls?)
    posPositions = detectBallBoundaries(roiBinarized, roiRect);
    
    if ~isempty(posPositions)
        if length(posPositions(:,1)) == 1
            ball_found = true;
            ball_position = posPositions;
        end
    end
end

if ball_found == true
    positions = [positions, ball_position];
end

% here we start displaying the frames..
if(ball_found_prev)
    drawTrackingWithPrevious(frame, positions);
else
    drawTrackingWithoutPrevious(frame, positions, roiRect);
    
end

frame_previous = frame;
ball_found_prev = ball_found;

function centroid = detectBallBoundaries(roiBinarized, roiRect)
    [B, L] = bwboundaries(roiBinarized);
    %posPositions = [];
    %figure; imshow(frame);
    %hold on;
    centroid = [];
    
    stats = regionprops(L,'Area','Centroid');

    currentMax = 0.7; % let's see how to set this threshold better 

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
        if metric > currentMax
            centroid = stats(k).Centroid;
            centroid(1) = centroid(1) + roiRect(1);
            centroid(2) = centroid(2) + roiRect(2);
            %posPositions = [posPositions, centroid];
            currentMax = metric;
            %plot(centroid(1),centroid(2),'ko');
        end

      %text(boundary(1,2)-35,boundary(1,1)+13,metric_string,'Color','y',...
       %    'FontSize',14,'FontWeight','bold');

    end

end

function roi = calcRoiSize(position, frame_size, ROI_WIDTH, ROI_HEIGHT)
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
    
    roi = [x_min, y_min, x_max - x_min, y_max - y_min ];    
end

function roi_binarized = thresholdImg(inputRoi, BALL_COLOR_HSV_MIN, BALL_COLOR_HSV_MAX, CLOSING_KERNEL_LENGTH, SMALLESTACCEPTABLEAREA)
    % REMEBMER: i wanna check the possibility of not calculating the edges and always use - 255
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
    % Keep areas only if they're bigger than this. THis needs to be tuned on the fly.

    
    ballColorMask = uint8(bwareaopen(ballColorMask, SMALLESTACCEPTABLEAREA));
    
    % morphological closing operation, imclose().
    structuringElement = strel('disk', CLOSING_KERNEL_LENGTH); % TAKE A LOOK HERE!
    ballColorMask = imclose(ballColorMask, structuringElement);
    
    % Fill in any holes
    ballColorMask = uint8(imfill(ballColorMask, 'holes'));
    
    % convert the type of ballColorMask to the same data type as hueBand.
    roi_binarized = cast(ballColorMask, class(hueBand));
    
    % Use the ball color object mask to mask out the ball color-only portions of the rgb image.
    %maskedImageH = ballColorMask .* hueBand;
    %maskedImageS = ballColorMask .* saturationBand;
    %maskedImageI = ballColorMask .* intensityBand;
    
    % Concatenate the masked color bands to form the rgb image.
    %roi_binarized = cat(3, maskedImageH, maskedImageS, maskedImageI);
    
    %imshow(roi_binarized);
end

function [circles, radii] = getCandidateBalls(inputRoi, roiRect, RADIUSRANGE)
    cannyRoi = edge(inputRoi, 'Canny');
    
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
        
    %figure;
    %imshow(frame);
    %hold on;
    
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

