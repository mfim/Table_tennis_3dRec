clc;
clear;
close all;
fontSize = 20;

MAX_ITERATIONS = 3;
BALL_SIZE = 22; % pixels of the ball 
% just for now.. keep the ballcolor
ballColor = 'o';
first_threshold = 3; second_threshold = 8;

%firstPosition = 'Video/5-m-60fps.mp4';
firstPosition = 'Video/5-m-120fps.mp4';
%secondPosition = 'Video/5-r-60fps.mp4';
secondPosition = 'Video/5-r-120fps.mp4';
%videoName = 'Video/1-m-30fps.mp4';


% prompt = 'Frames to be skipped: ';
% skips = str2double(input(prompt, 's'));
% if isnan(skips) || fix(skips) ~= skips
%   disp('Please enter an integer');
%   return;
% end

intrinsics = Hero3_120fps();
vFirst = VideoReader(firstPosition, 'CurrentTime', 0.3);
vSecond = VideoReader(secondPosition, 'CurrentTime', 0.3);
frameFirstCam = readFrame(vFirst);
frameSecondCam = readFrame(vSecond);
[paramsFirstCam, rotationMatrixFirstCam, translationVectorFirstCam] = calcExtrinsics(frameFirstCam, intrinsics);      
[paramsSecondCam, rotationMatrixSecondCam, translationVectorSecondCam] = calcExtrinsics(frameSecondCam, intrinsics);
clear vFirst vSecond;

[PositionsFirstCam, CurveFirstCam, BounceTsFirstCam, BounceCoordFirstCam, StrikeTsFirstCam, StrikeCoordFirstCam] = ballTracking(intrinsics, firstPosition, ... 
     MAX_ITERATIONS, BALL_SIZE, ballColor, first_threshold, second_threshold, 30);
[PositionsSecondCam, CurveSecondCam, BounceTsSecondCam, BounceCoordSecondCam, StrikeTsSecondCam, StrikeCoordSecondCam] = ballTracking(intrinsics, secondPosition,...
     MAX_ITERATIONS, BALL_SIZE, ballColor, first_threshold, second_threshold, 0);

 % not working currently
offset = syncCam(BounceCoordFirstCam, BounceCoordSecondCam, BounceTsFirstCam, BounceTsSecondCam, ...
    paramsFirstCam, rotationMatrixFirstCam, translationVectorFirstCam, ...
    paramsSecondCam, rotationMatrixSecondCam, translationVectorSecondCam);

% zeroTime = CurveFirstCam(1,3);
% CurveFirstCam(:,3) = CurveFirstCam(:,3) - zeroTime; 
CurveSecondCam(:,3) = CurveSecondCam(:,3) - (offset);
%CurveSecondCam = CurveSecondCam(:,3) > 0;

% Handle the start 
if(CurveFirstCam(1,3) > CurveSecondCam(1,3))
    CurveSecondCam = CurveSecondCam(CurveSecondCam(:,3) >= CurveFirstCam(1,3), :);
else
    CurveFirstCam = CurveFirstCam(CurveFirstCam(:,3) >= CurveSecondCam(1,3), :);
end

% Handle the end 
if(CurveFirstCam(length(CurveFirstCam(:,1)) ,3) < CurveSecondCam(length(CurveFirstCam(:,1)),3))
    CurveSecondCam = CurveSecondCam(CurveSecondCam(:,3) < CurveFirstCam(length(CurveFirstCam(:,1)),3), :);
else
    CurveFirstCam = CurveFirstCam(CurveFirstCam(:,3) < CurveSecondCam(length(CurveFirstCam(:,1)),3),:);
end

MatchSecondCam = interp1(CurveSecondCam(:, 3), CurveSecondCam(:, 1:2), CurveFirstCam(:,3) , 'spline');

Trajectory = triangulate(CurveFirstCam(:, 1:2), MatchSecondCam, paramsFirstCam, paramsSecondCam);
%Trajectory = triangulate(CurveFirstCam(1:225, 1:2), CurveSecondCam(:, 1:2), paramsFirstCam, paramsSecondCam);


plotTrajectory(Trajectory);


function plotTrajectory(Trajectory)
figure;

curve = animatedline('lineWidth', 2);
set(gca, 'XLim', [-50 280], 'YLim', [-80 280], 'ZLim', [-50 100]);
view(43, 24);
plotTable(min(abs(Trajectory(:,3))));
grid on;
hold on; 

for i=1:length(Trajectory)
  addpoints(curve, Trajectory(i, 1), Trajectory(i, 2), -Trajectory(i, 3));
  head = scatter3(Trajectory(i, 1), Trajectory(i, 2), -Trajectory(i, 3),'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'y');
  drawnow
  pause(0.01);
  delete(head);
end

end


%Camera Calibration Matrix - Intrinsics:
% GOPRO Hero3 720p - 120fps - Narrow
function [intrinsics] = Hero3_120fps()

focalLength = [1101, 1101];
principalPoint = [639.5, 359.5];
imageSize = [720,1280];
radialDistortion = [-0.359, 0.279];

intrinsics = cameraIntrinsics(focalLength, principalPoint, imageSize, 'RadialDistortion', radialDistortion);

%j = undistortImage(frame, intrinsics);
%j = undistortImage(frame, intrinsics, 'OutputView', 'full');
%figure; imshowpair(frame,j,'montage');

end

function [cameraParameters, rotationMatrix, translationVector] = calcExtrinsics(frame, intrinsics)

% official table dimension in cm
worldPoints = [0, 0; 152.5, 0; 0, 274; 152.5, 274];

% rectify image
j = undistortImage(frame, intrinsics, 'OutputView', 'full');
figure; imshow(j, []);

% input the corners on rectified image
% ALWAYS A Z
[x, y] = getpts;
imagePoints = [x, y];

[rotationMatrix, translationVector] = extrinsics(imagePoints, worldPoints, intrinsics);
cameraParameters = cameraMatrix(intrinsics, rotationMatrix, translationVector);
%cameraParameters = cameraParameters('IntrinsicMatrix', intrinsics.IntrinsicMatrix, ...);

%tablePlane = BounceCoordSecondCam/paramsSecondCam;


% plot camera position 
[orientation, location] = extrinsicsToCameraPose(rotationMatrix, translationVector);

figure; plotCamera('Location', location, 'Orientation', orientation, 'Size', 20);
hold on;
pcshow([worldPoints, zeros(size(worldPoints, 1),1)], ...
    'VerticalAxisDir','down','MarkerSize',40);

end