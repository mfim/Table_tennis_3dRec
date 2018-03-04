% REMEBMER: i wanna check the possibility of not calculating the edges and always use - 255
inputRoi = imread('yellow_ball.jpg');

roi_hsv_aux = rgb2hsv(inputRoi);
roi_hsv = uint8(255* roi_hsv_aux);

hueBand = roi_hsv (:,:,1);
saturationBand = roi_hsv(:,:,2);
% this will be disconsidered
intensityBand = roi_hsv(:,:,3);

%debugging purposes
subplot(3, 4, 2);
imshow(hueBand);
title('hue Band', 'FontSize', fontSize);
subplot(3, 4, 3);
imshow(saturationBand);
title('saturation Band', 'FontSize', fontSize);
subplot(3, 4, 4);
imshow(intensityBand);
title('intensity Band', 'FontSize', fontSize);
message = sprintf('These are the individual HSV bands.\nNow we will compute the image histograms.');

% Compute and plot the hue historgram
hH = subplot(3, 4, 6);
[countsH, grayLevelsH] = imhist(hueBand);
maxGLValueH = find(countsH > 0, 1, 'last');
maxCountH = max(countsH);
bar(countsH, 'r');
grid on;
xlabel('Gray Levels');
ylabel('Pixel Count');
title('Histogram of Hue Band', 'FontSize', fontSize);

% Compute and plot the saturation historgram
hS = subplot(3, 4, 7);
[countsS, grayLevelsS] = imhist(saturationBand);
maxGLValueS = find(countsS > 0, 1, 'last');
maxCountS = max(countsS);
bar(countsS, 'g', 'BarWidth', 0.95);
grid on;
xlabel('Gray Levels');
ylabel('Pixel Count');
title('Histogram of Saturation Band', 'FontSize', fontSize);

% Compute and plot the intensity histogram
hI = subplot(3, 4, 8);
[countsI, grayLevelsI] = imhist(intensityBand);
maxGLValueI = find(countsI > 0, 1, 'last');
maxCountI = max(countsI);
bar(countsI, 'b');
grid on;
xlabel('Gray Levels');
ylabel('Pixel Count');
title('Histogram of intensity Band', 'FontSize', fontSize);

%set axes to same widhth and height
maxGL = max([maxGLValueH, maxGLValueS, maxGLValueI]);

maxCount = max([maxCountH,  maxCountS, maxCountI]);
axis([hH hS hI], [0 maxGL 0 maxCount]);

% Plot all 3 histograms in one plot.
subplot(3, 4, 5);
plot(grayLevelsH, countsH, 'r', 'LineWidth', 2);
grid on;
xlabel('Gray Levels');
ylabel('Pixel Count');
hold on;
plot(grayLevelsS, countsS, 'g', 'LineWidth', 2);
plot(grayLevelsI, countsI, 'b', 'LineWidth', 2);
title('Histogram of All Bands', 'FontSize', fontSize);
maxGrayLevel = max([maxGLValueH, maxGLValueS, maxGLValueI]);

% trim the axis to the max gray level on the bright end
xlim([0 maxGrayLevel]);

%PlaceThresholdBars(6, BALL_COLOR_HSV_MIN(1), BALL_COLOR_HSV_MAX(1));
%PlaceThresholdBars(7, BALL_COLOR_HSV_MIN(2), BALL_COLOR_HSV_MAX(2));
%PlaceThresholdBars(8, BALL_COLOR_HSV_MIN(3), BALL_COLOR_HSV_MAX(3));

% Now apply each color band's particular thresholds to the color band
hueMask = (hueBand >= BALL_COLOR_HSV_MIN(1)) & (hueBand <= BALL_COLOR_HSV_MAX(1));
saturationMask = (saturationBand >= BALL_COLOR_HSV_MIN(2)) & (saturationBand <= BALL_COLOR_HSV_MAX(2));
intensityMask = (intensityBand >= BALL_COLOR_HSV_MIN(3)) & (intensityBand <= BALL_COLOR_HSV_MAX(3));

% Display the thresholded binary images.
fontSize = 16;
subplot(3, 4, 10);
imshow(hueMask, []);
title('Is-Hue Mask', 'FontSize', fontSize);
subplot(3, 4, 11);
imshow(saturationMask, []);
title('Is-Not-Saturation Mask', 'FontSize', fontSize);
subplot(3, 4, 12);
imshow(intensityMask, []);
title('Is-Not-intensity Mask', 'FontSize', fontSize);
% Combine the masks to find where all 3 are "true."
% Then we will have the mask of only the red parts of the image.
ballColorMask = uint8(hueMask & intensityMask & saturationMask);
subplot(3, 4, 9);
imshow(ballColorMask, []);
caption = sprintf('Mask of Only\nThe Ball colored Objects');
title(caption, 'FontSize', fontSize);

% let's filter small objects 
smallestAcceptableArea = 5000; % Keep areas only if they're bigger than this. THis needs to be tuned on the fly. 
% the ball in the test pict is REALLY big. remember to change this.
   
figure;
% Maximize the figure.
set(gcf, 'Position', get(0, 'ScreenSize'));

ballColorMask = uint8(bwareaopen(ballColorMask, smallestAcceptableArea));
subplot(3, 3, 1);
imshow(ballColorMask, []);
fontSize = 13;
caption = sprintf('bwareaopen() removed objects\nsmaller than %d pixels', smallestAcceptableArea);
title(caption, 'FontSize', fontSize);

% morphological closing operation, imclose().
structuringElement = strel('disk', CLOSING_KERNEL_LENGTH); % TAKE A LOOK HERE!
ballColorMask = imclose(ballColorMask, structuringElement);
subplot(3, 3, 2);
imshow(ballColorMask, []);
fontSize = 16;
title('Border smoothed', 'FontSize', fontSize);

% Fill in any holes
ballColorMask = uint8(imfill(ballColorMask, 'holes'));
subplot(3, 3, 3);
imshow(ballColorMask, []);
title('Regions Filled', 'FontSize', fontSize);

% convert the type of ballColorMask to the same data type as hueBand.
ballColorMask = cast(ballColorMask, class(hueBand));


% Use the red object mask to mask out the red-only portions of the rgb image.
maskedImageH = ballColorMask .* hueBand;
maskedImageS = ballColorMask .* saturationBand;
maskedImageI = ballColorMask .* intensityBand;
% Show the masked off hue image.
subplot(3, 3, 4);
imshow(maskedImageH);
title('Masked Hue Image', 'FontSize', fontSize);
% Show the masked off saturation image.
subplot(3, 3, 5);
imshow(maskedImageS);
title('Masked Saturation Image', 'FontSize', fontSize);
% Show the masked off intensity image.
subplot(3, 3, 6);
imshow(maskedImageI);
title('Masked intensity Image', 'FontSize', fontSize);
% Concatenate the masked color bands to form the rgb image.
maskedRGBImage = cat(3, maskedImageH, maskedImageS, maskedImageI);
% Show the masked off, original image.
subplot(3, 3, 8);
imshow(maskedRGBImage);
fontSize = 13;
caption = sprintf('Masked Original Image\nShowing Only the ball color Objects');
title(caption, 'FontSize', fontSize);
% Show the original image next to it.
subplot(3, 3, 7);
imshow(inputRoi);
title('The Original Image (Again)', 'FontSize', fontSize);