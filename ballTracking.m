function [final_positions, vq, bounceTs, bouncePoint, strikeTs, strikePoint] = ballTracking(intrinsics, videoName, MAX_ITERATIONS, BALL_SIZE, ballColor, first_threshold, second_threshold, skips)
% ---- BALL RELATED CONSTANTS ---- %
BALL_AREA = pi * (BALL_SIZE/2)^2;
BALL_PERIMETER = 2 * pi * (BALL_SIZE/2);
RADIUSRANGE = [round(BALL_SIZE/2-0.3*BALL_SIZE/2), round(BALL_SIZE/2+0.3*BALL_SIZE/2)];
ROI_WIDTH = BALL_SIZE * 10; %this can vary with the video FPS
ROI_HEIGHT = BALL_SIZE * 10;
ADAPTIVE_WIDTH = ROI_WIDTH;
ADAPTIVE_HEIGHT = ROI_HEIGHT;

% kernel for blurring. VALUE: 3-sigma rule
BLUR_KERNEL_LENGTH = (BALL_SIZE-1)/6;
CLOSING_KERNEL_LENGTH = round(BALL_SIZE/2);

% ---- END OF BALL RELATED CONSTANTS ---- %

% --- var initialization --- %
skipped = 0;
%SLOPES = [];
positions = [];
ball_found_prev = false;

first_frame = true;
% to debug iterations
debug = 0;

v = VideoReader(videoName);

% rough estimation of time in frames
% VideoReader does have a methdod currentTime, but
% for some reason it was giving wron information
numFrames = round(v.duration * v.FrameRate);
timeSingleFrame = v.duration / numFrames;
currentTime = 0;


% --- MAIN LOOP --- %
while hasFrame(v)
    
    % skip the bad frames at the beginning
    if skipped < skips
        if hasFrame(v)
            frame = readFrame(v);
            frame = undistortImage(frame, intrinsics);
        end
        skipped=skipped+1;
    else
        frame = readFrame(v);
        frame = undistortImage(frame, intrinsics);
        ball_found = 0;
        candidateRoiRect = [];
        candidateBallsRoi = [];
        
        % calculate the ROI
        if ball_found_prev == 1
            ADAPTIVE_WIDTH = ROI_WIDTH;
            ADAPTIVE_HEIGHT = ROI_HEIGHT;
            roiRect = calcRoiSize(positions(length(positions(:,1)), :), [size(frame,2), size(frame,1)], ROI_WIDTH, ROI_HEIGHT);
            roi = imcrop(frame, roiRect);
        else
            if first_frame  == true
                figure;
                imshow(frame);
                [positions(1), positions(2)] = getpts;
                positions(3) = currentTime;
                [FIRST_PASS_HSV_MAX, FIRST_PASS_HSV_MIN, SECOND_PASS_HSV_MAX, SECOND_PASS_HSV_MIN, h] = hsvRangesDef(positions(1), positions(2), ...
                    frame, first_threshold, second_threshold);
                roiRect = calcRoiSize(positions, [size(frame,2), size(frame,1)], ROI_WIDTH, ROI_HEIGHT);
                frame = readFrame(v);
                frame = undistortImage(frame, intrinsics);
                currentTime = currentTime + timeSingleFrame;
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
        
        debug = debug+1;
        if (debug == 190)
            disp('dummy');
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
                candidateRoiRect(i,:) = calcRoiSize(stats(i).Centroid, [size(frame,2), size(frame,1)], 2.5*BALL_SIZE, 2.5*BALL_SIZE);
                candidateBallsRoi{i} = imcrop(roiBlurred, candidateRoiRect(i, :));
            end
        else
            candidateRoiRect = roiRect;
            candidateBallsRoi{1} = roiBlurred;
        end
        
        posError = [0, 0, inf];
        
        if ~isempty(stats)
            
            % check the returns here, seem wrong!
            % --- SECOND PASS: RELAXED THRESHOLD --- %
            [candidateBallsRoiBin, candidateBallsHSV] = secondPass(candidateBallsRoi, ...
                ballColor, h, MAX_ITERATIONS, second_threshold, SECOND_PASS_HSV_MIN, SECOND_PASS_HSV_MAX, CLOSING_KERNEL_LENGTH, ...
                BALL_AREA);
            
            
            % --- FIND THE BALL! --- %
            for i = 1:length(candidateBallsRoiBin)
                candidateRoiBin = candidateBallsRoiBin{i};
                
                % in case we're using a ROI inside the ROI
                if ~isequal(candidateRoiRect, roiRect)
                    candidateRoiRect(i, 1) = candidateRoiRect(i, 1) + roiRect(1);
                    candidateRoiRect(i, 2) = candidateRoiRect(i, 2) + roiRect(2);
                end
                
                %[error, SLOPES ]= twoStageBallDetection(candidateRoiBin, candidateRoiRect(i,:), positions, currentTime, BALL_SIZE, BALL_AREA, BALL_PERIMETER, SLOPES);
                error = twoStageBallDetection(candidateRoiBin, candidateRoiRect(i,:), positions, currentTime, BALL_SIZE, BALL_AREA, BALL_PERIMETER);
                
                if error < posError(3)
                    posError = error;
                    SECOND_PASS_HSV_MIN = candidateBallsHSV(i, 1:3);
                    SECOND_PASS_HSV_MAX = candidateBallsHSV(i, 4:6);
                    second_threshold = candidateBallsHSV(i, 7);
                end
                
            end
            
            if posError(3) < inf
                ball_found = true;
                positions = [positions ; posError(1), posError(2), currentTime];
            end
        end
        frame_previous = frame;
        ball_found_prev = ball_found;
    end
    
    currentTime = currentTime + timeSingleFrame;
end
% ---- END OF MAIN LOOP ---- %
[~, outlier] = hampel(positions, 1);
final_positions = [];
for i = 1:length(outlier)
    if ~(outlier(i,1) || outlier(i,2) || outlier(i,3) == 1)
        final_positions = [final_positions; positions(i,:)];
    end
end

% let's plot the balls found, include the radius
hold on; plot(final_positions(:,1), final_positions(:,2), 'ro');
% Change to 0:001 if video is 120fps
interval = final_positions(1,3):0.01:final_positions(length(final_positions),3);
vq = interp1(final_positions(:, 3), final_positions(:, 1:2), interval , 'spline');
vq(:, 3) = interval;
hold on; plot(vq(:,1), vq(:,2), 'g--');

% search for the bounce
% analyzing y in 5 frames
candidateTs = [];
for i= 4:length(final_positions)-3
    if ((final_positions(i-3, 2) < final_positions(i-2, 2)) && ...
        (final_positions(i-2, 2) < final_positions(i-1, 2)) && (final_positions(i-1, 2) < final_positions(i,2)) && ...
            (final_positions(i, 2) > final_positions(i+1, 2)) && (final_positions(i+1, 2) > final_positions(i+2,2)) ...
            && (final_positions(i+2, 2) > final_positions(i+3,2)))
        % in a bounce there is no change in the x direction 
            if((final_positions(i-2, 1) > final_positions(i-1, 1)) && final_positions(i-1, 1) > final_positions(i,1) && ... 
                    final_positions(i,1) > final_positions(i+1,1) && final_positions(i+1, 1) > final_positions(i+2,1) || ...
                    (final_positions(i-2, 1) < final_positions(i-1, 1)) && final_positions(i-1, 1) < final_positions(i,1) && ... 
                    final_positions(i,1) < final_positions(i+1,1) && final_positions(i+1, 1) < final_positions(i+2,1)) 
                candidateTs = [candidateTs; i];
            end
    end
end

% estimate the inside the frames point of bounce
bounceTs = [];
bouncePoint = [];
for j = 1:length(candidateTs)
    interval = final_positions(candidateTs(j)-2,3):0.001:final_positions(candidateTs(j)+2,3);
    vq2 = interp1(final_positions(:, 3), final_positions(:, 1:2), interval , 'spline');
    for i = 4:length(vq2)-3
        if ((vq2(i-3, 2) < vq2(i-2, 2)) && (vq2(i-2, 2) < vq2(i-1, 2)) && (vq2(i-1, 2) < vq2(i,2)) && ...
                (vq2(i, 2) > vq2(i+1, 2)) && (vq2(i+1, 2) > vq2(i+2,2)) && (vq2(i+2, 2) > vq2(i+3,2)))
            bounceTs = [bounceTs ; interval(i)];
            bouncePoint = [bouncePoint; vq2(i, :)];
            
        end
    end
end
hold on; plot(bouncePoint(:,1), bouncePoint(:,2), 'w*');


% search for the strike
% analyzing y in 5 frames
candidateTs = [];
for i= 3:length(final_positions)-2
    if ((final_positions(i-2, 1) < final_positions(i-1, 1)) && (final_positions(i-1, 1) < final_positions(i,1)) && ...
            (final_positions(i, 1) > final_positions(i+1, 1)) && (final_positions(i+1, 1) > final_positions(i+2,1)) || ...
            (final_positions(i-2, 1) > final_positions(i-1, 1)) && (final_positions(i-1, 1) > final_positions(i,1)) && ...
            (final_positions(i, 1) < final_positions(i+1, 1)) && (final_positions(i+1, 1) < final_positions(i+2,1)))
        candidateTs = [candidateTs; i];
    end
end

% estimate the inside the frames point of strike
strikeTs = [];
strikePoint = [];
for j = 1:length(candidateTs)
    interval = final_positions(candidateTs(j)-2,3):0.001:final_positions(candidateTs(j)+2,3);
    vq2 = interp1(final_positions(:, 3), final_positions(:, 1:2), interval , 'spline');
    for i = 4:length(vq2)-3
        if (((vq2(i-3, 1) < vq2(i-2, 1)) && (vq2(i-2, 1) < vq2(i-1, 1)) && (vq2(i-1, 1) < vq2(i,1)) && ...
                (vq2(i, 1) > vq2(i+1, 1)) && (vq2(i+1, 1) > vq2(i+2,1)) && (vq2(i+2, 1) > vq2(i+3,1))) || ...
            (vq2(i-3, 1) > vq2(i-2, 1)) && (vq2(i-2, 1) > vq2(i-1, 1)) && (vq2(i-1, 1) > vq2(i,1)) && ...
                (vq2(i, 1) < vq2(i+1, 1)) && (vq2(i+1, 1) < vq2(i+2,1)) && (vq2(i+2, 1) < vq2(i+3,1)))
            strikeTs = [strikeTs ; interval(i)];
            strikePoint = [strikePoint; vq2(i, :)];
            
        end
    end
end

end


function [FIRST_PASS_HSV_MAX, FIRST_PASS_HSV_MIN, SECOND_PASS_HSV_MAX, SECOND_PASS_HSV_MIN, h] = hsvRangesDef(x, y, frame, first_threshold, second_threshold)

% for a better accuracy we average the color with its neighbour pixels
roiBallRect = calcRoiSize([x y],[size(frame,2), size(frame,1)], 5, 5);
roiBall = imcrop(frame, roiBallRect);
% Separate to RGB channel
Ir = roiBall(:,:,1);
Ig = roiBall(:,:,2);
Ib = roiBall(:,:,3);

% Calculate average RGB of the region
Rave = round(double(mean(Ir(:))));
Gave = round(double(mean(Ig(:))));
Bave = round(double(mean(Ib(:))));

ball_hsv = rgb2hsv([Rave/255,Gave/255, Bave/255]);
h = ball_hsv(:,1)*360;


FIRST_PASS_HSV_MIN = [ball_hsv(:,1) - first_threshold/360, 45/100, 70/100];
FIRST_PASS_HSV_MAX = [ball_hsv(:,1) + first_threshold/360, 1, 1];
% may be add a variable in the saturation as well
SECOND_PASS_HSV_MIN = [ball_hsv(:,1) - second_threshold/360, 45/100, 70/100];
SECOND_PASS_HSV_MAX = [ball_hsv(:,1) + second_threshold/360, 1, 1];

% handle edges
if(FIRST_PASS_HSV_MIN(1) < 0)
    FIRST_PASS_HSV_MIN(1) = 0;
end
if(FIRST_PASS_HSV_MAX(1) > 1)
    FIRST_PASS_HSV_MAX(1) = 1;
end
if(SECOND_PASS_HSV_MIN(1) < 0)
    SECOND_PASS_HSV_MIN(1) = 0;
end
if(SECOND_PASS_HSV_MAX(1) > 1)
    SECOND_PASS_HSV_MAX(1) = 1;
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

function [xc, yc, R] = CircleFitByPratt(x,y)

%--------------------------------------------------------------------------
%
%     Circle fit by Pratt
%      V. Pratt, "Direct least-squares fitting of algebraic surfaces",
%      Computer Graphics, Vol. 21, pages 145-152 (1987)
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: Par = [a b R] is the fitting circle:
%                           center (a,b) and radius R
%
%     Note: this fit does not use built-in matrix functions (except "mean"),
%           so it can be easily programmed in any programming language
%
%--------------------------------------------------------------------------

n = size(x,1);      % number of data points

XY(:,1) = x;
XY(:,2) = y;

centroid = mean(XY);   % the centroid of the data set

%     computing moments (note: all moments will be normed, i.e. divided by n)

Mxx=0; Myy=0; Mxy=0; Mxz=0; Myz=0; Mzz=0;

for i=1:n
    Xi = XY(i,1) - centroid(1);  %  centering data
    Yi = XY(i,2) - centroid(2);  %  centering data
    Zi = Xi*Xi + Yi*Yi;
    Mxy = Mxy + Xi*Yi;
    Mxx = Mxx + Xi*Xi;
    Myy = Myy + Yi*Yi;
    Mxz = Mxz + Xi*Zi;
    Myz = Myz + Yi*Zi;
    Mzz = Mzz + Zi*Zi;
end

Mxx = Mxx/n;
Myy = Myy/n;
Mxy = Mxy/n;
Mxz = Mxz/n;
Myz = Myz/n;
Mzz = Mzz/n;

%    computing the coefficients of the characteristic polynomial

Mz = Mxx + Myy;
Cov_xy = Mxx*Myy - Mxy*Mxy;
Mxz2 = Mxz*Mxz;
Myz2 = Myz*Myz;

A2 = 4*Cov_xy - 3*Mz*Mz - Mzz;
A1 = Mzz*Mz + 4*Cov_xy*Mz - Mxz2 - Myz2 - Mz*Mz*Mz;
A0 = Mxz2*Myy + Myz2*Mxx - Mzz*Cov_xy - 2*Mxz*Myz*Mxy + Mz*Mz*Cov_xy;
A22 = A2 + A2;

epsilon=1e-12;
ynew=1e+20;
IterMax=20;
xnew = 0;

%    Newton's method starting at x=0

for iter=1:IterMax
    yold = ynew;
    ynew = A0 + xnew*(A1 + xnew*(A2 + 4.*xnew*xnew));
    if (abs(ynew)>abs(yold))
        disp('Newton-Pratt goes wrong direction: |ynew| > |yold|');
        xnew = 0;
        break;
    end
    Dy = A1 + xnew*(A22 + 16*xnew*xnew);
    xold = xnew;
    xnew = xold - ynew/Dy;
    if (abs((xnew-xold)/xnew) < epsilon), break, end
    if (iter >= IterMax)
        disp('Newton-Pratt will not converge');
        xnew = 0;
    end
    if (xnew<0.)
        fprintf(1,'Newton-Pratt negative root:  x=%f\n',xnew);
        xnew = 0;
    end
end

%    computing the circle parameters

DET = xnew*xnew - xnew*Mz + Cov_xy;
Center = [Mxz*(Myy-xnew)-Myz*Mxy , Myz*(Mxx-xnew)-Mxz*Mxy]/DET/2;

xc = Center(:,1) +centroid(:,1);
yc = Center(:,2) +centroid(:,2);
R = sqrt(Center*Center'+Mz+2*xnew);
%Par = [Center+centroid , sqrt(Center*Center'+Mz+2*xnew)];

end    %    CircleFitByPratt


function minError= twoStageBallDetection(roiBinarized, roiRect, positions, nextTs, BALL_SIZE, BALL_AREA, BALL_PERIMETER)
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
        [xc, yc, R] = CircleFitByPratt(upperContour(:,2), upperContour(:,1));
        
        % calculate the error
        n = length(upperContour(:,1));
        sumError = 0;
        for i=1:n
            d = sqrt((upperContour(i,2)- xc)^2 + (upperContour(i,1) - yc)^2);
            sumError = sumError + abs((d - R)/R);
        end
        
        Eruc = sumError/n;
        % fine tune the threshold
        if (Eruc < 0.1 && R < 0.8*BALL_SIZE && (round(yc) > 0 && round(yc) < length(roiBinarized(:,1)) && ...
                round(xc) > 0 && round(xc) < length(roiBinarized(:,2)) && roiBinarized(round(yc), round(xc))))
            RUC = true;
        else
            
            lowerContour = boundary( boundary(:,1) >= stats(k).BoundingBox(4)/2 + stats(k).BoundingBox(2), :);
            % fit in a circle
            if length(upperContour) < 4
                RUC = false;
            else
                [xc, yc, R] = CircleFitByPratt(lowerContour(:,2), lowerContour(:,1));
                
                % calculate the error
                n = length(lowerContour(:,1));
                sumError = 0;
                for i=1:n
                    d = sqrt((lowerContour(i,2)- xc)^2 + (lowerContour(i,1) - yc)^2);
                    sumError = sumError + abs((d - R)/R);
                end
                
                Eruc = sumError/n;
                
                if (Eruc < 0.1 && R < 0.8*BALL_SIZE && (round(yc) > 0 && round(yc) < length(roiBinarized(:,1)) && ...
                        round(xc) > 0 && round(xc) < length(roiBinarized(:,2))&& roiBinarized(round(yc), round(xc))))
                    RUC = true;
                else
                    RUC = false;
                end
            end
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
            % the 0.2 could vary with the frameRate
            if (distC > 0.1 * BALL_SIZE && distC < 12 * BALL_SIZE) % set this threshold!
                if((positions(length(positions(:,1)), 1) < positions(length(positions(:,1))-1, 1) && ...
                        centroid(:,1) > positions(length(positions(:,1)), 1)) || ...
                        positions(length(positions(:,1)), 1) > positions(length(positions(:,1))-1, 1) && ...
                        centroid(:,1) < positions(length(positions(:,1)), 1))
                    if distC < 5 * BALL_SIZE
                        mc = true;
                    end
                else
                    mc = true;
                end
                mc = true;
            end
            
            %                 distP = pdist([prediction; positions(length(positions(:,1)), 1:2)], 'euclidean');
            %                 if distP > 0.5 * BALL_SIZE
            %                     mp = true;
            %                 end
        elseif length(positions(:,1)) == 1
            % still don't have enough for a prediction, check only
            % movement
            distC = pdist([centroid; positions(1,1:2)], 'euclidean');
            if distC > 0.2 * BALL_SIZE % set this threshold!
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
        if (first_threshold > 2.5)
            first_threshold = first_threshold*0.90;
            
            if strcmp(ballColor, 'o')
                FIRST_PASS_HSV_MIN(1) = (h - first_threshold)/360;
                %FIRST_PASS_HSV_MIN(2) = (75 - first_threshold/2)/100;
                FIRST_PASS_HSV_MAX(1) = (h + first_threshold)/360;
            else
                %FIRST_PASS_HSV_MIN = [0, 0, 0.6];
                FIRST_PASS_HSV_MAX(2) = h+first_threshold/100;
            end
        end
    elseif candidateBallsInfo.NumObjects == 1
        candidate_found = true;
        best_first_pass = [iteration,candidateBallsInfo.NumObjects];
        BEST_FIRST_PASS_HSV_MIN =  FIRST_PASS_HSV_MIN;
        BEST_FIRST_PASS_HSV_MAX =  FIRST_PASS_HSV_MAX;
        best_threshold = first_threshold;
    else
        % decrease the threshold
        if (first_threshold < 20)
            first_threshold = first_threshold * 1.1;
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
                %                 else
                %                     %FIRST_PASS_HSV_MIN = [0, 0, 0.60];
                %                     FIRST_PASS_HSV_MAX(2) = h+first_threshold/100;
                %                     if(FIRST_PASS_HSV_MAX(2) > 1)
                %                         FIRST_PASS_HSV_MAX(2) = 1;
                %                     end
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
            if  aux_area < 0.85
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
                if (it_threshold > 5)
                    it_threshold = it_threshold*0.90;
                    
                    if strcmp(ballColor, 'o')
                        IT_HSV_MIN(1) = (h - it_threshold)/360;
                        %SECOND_PASS_HSV_MIN(2) = (55 - second_threshold/2)/100;
                        IT_HSV_MAX(1) = (h + it_threshold)/360;
                    else
                        IT_HSV_MAX(2) = h + it_threshold/100;
                    end
                end
            else
                
                if best_it(2) > BALL_AREA - area
                    best_it = [iteration, BALL_AREA - area];
                    BEST_ITERATION_HSV_MAX = IT_HSV_MAX;
                    BEST_ITERATION_HSV_MIN = IT_HSV_MIN;
                    best_iteration_threshold = it_threshold;
                end
                
                % decrease the threshold
                it_threshold = it_threshold * 1.1;
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
                %ends
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
    
    candidateBallsRoi{i} = candidateRoiBinarized;
    % EITHER KEEP BEST_IT OR THIS BELOW!!
    if best_it(2) == inf
        candidateBallsHSV = [candidateBallsHSV; SECOND_PASS_HSV_MIN, SECOND_PASS_HSV_MAX, second_threshold];
    else
        candidateBallsHSV = [candidateBallsHSV; BEST_ITERATION_HSV_MIN, BEST_ITERATION_HSV_MAX, best_iteration_threshold];
    end
end

end
