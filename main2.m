%%
clc;
clear;
close all;
fontSize = 20;

MAX_ITERATIONS = 3;
BALL_SIZE = 22; % pixels of the ball 
% just for now.. keep the ballcolor
ballColor = 'o';
first_threshold = 3; second_threshold = 8;

% File name definition 
firstPosition = 'Video/5-m-120fps.mp4';
secondPosition = 'Video/5-r-120fps.mp4';


% 
v1 = VideoReader(firstPosition, 'CurrentTime', 0.3); %1.5 for 120-m 3.8 for 120-r
v2 = VideoReader(secondPosition, 'CurrentTime', 0.8); %0.5 for 60-m-c 1 for 60-r-c
                                                    %3 for 60-m and 3 for
                                                    %60-r
% 
camera1 = cb_cab_r(121); %Getting the camera matrice given the fps 
camera2 = cb_cab_r(120);

% Get the frame to detect the table on 
view1 = readFrame(v1);
view2 = readFrame(v2);

% Remove the radial distortion before 
view1 = undistortImage(view1, camera1);
view2 = undistortImage(view2, camera2);

[~, corners1] = table_detection(view1 ,1);
[~, corners2] = table_detection(view2 ,1);
%%
% Need to reorder the points to match the real table coordinates
cornersr1 = [
    corners1(find(corners1(:,2) == max(corners1(:,2))),:);
    corners1(find(corners1(:,1) == min(corners1(:,1))),:);
    corners1(find(corners1(:,1) == max(corners1(:,1))),:);
    corners1(find(corners1(:,2) == min(corners1(:,2))),:);
];

cornersr2 = [
    corners2(find(corners2(:,1) == min(corners2(:,1))),:);
    corners2(find(corners2(:,2) == min(corners2(:,2))),:);
    corners2(find(corners2(:,2) == max(corners2(:,2))),:);
    corners2(find(corners2(:,1) == max(corners2(:,1))),:);
];


[paramsFirstCam, rotationMatrixFirstCam, translationVectorFirstCam] = calcExtrinsicsTemp(view1, camera1, cornersr1);
[paramsSecondCam, rotationMatrixSecondtCam, translationVectorSecondCamn] = calcExtrinsicsTemp(view2, camera2, cornersr2);
%[paramsFirstCam, rotationMatrixFirstCam, translationVectorFirstCam] = calcExtrinsics(view1, camera1);      
%[paramsSecondCam, rotationMatrixSecondCam, translationVectorSecondCam] = calcExtrinsics(view2, camera2);
clear v1, v2;


vFirst = VideoReader(firstPosition, 'CurrentTime', 0.3);
vSecond = VideoReader(secondPosition, 'CurrentTime', 0.3);


% track ball
[PositionsFirstCam, CurveFirstCam, BounceTsFirstCam, BounceCoordFirstCam, StrikeTsFirstCam, StrikeCoordFirstCam] = ballTracking(camera1, firstPosition, ... 
     MAX_ITERATIONS, BALL_SIZE, ballColor, first_threshold, second_threshold, 30);
[PositionsSecondCam, CurveSecondCam, BounceTsSecondCam, BounceCoordSecondCam, StrikeTsSecondCam, StrikeCoordSecondCam] = ballTracking(camera2, secondPosition,...
     MAX_ITERATIONS, BALL_SIZE, ballColor, first_threshold, second_threshold, 0);

offset = syncCam(BounceCoordFirstCam, BounceCoordSecondCam, BounceTsFirstCam, BounceTsSecondCam, paramsFirstCam, paramsSecondCam);

CurveSecondCam(:,3) = CurveSecondCam(:,3) - (offset);

% Handle the start 
if(CurveFirstCam(1,3) > CurveSecondCam(1,3))
    CurveSecondCam = CurveSecondCam(CurveSecondCam(:,3) >= CurveFirstCam(1,3), :);
else
    CurveFirstCam = CurveFirstCam(CurveFirstCam(:,3) >= CurveSecondCam(1,3), :);
end

% Handle the end 
if(CurveFirstCam(length(CurveFirstCam(:,1)) ,3) < CurveSecondCam(length(CurveSecondCam(:,1)),3))
    CurveSecondCam = CurveSecondCam(CurveSecondCam(:,3) < CurveFirstCam(length(CurveFirstCam(:,1)),3), :);
else
    CurveFirstCam = CurveFirstCam(CurveFirstCam(:,3) < CurveSecondCam(length(CurveSecondCam(:,1)),3),:);
end

% estimate same points
MatchSecondCam = interp1(CurveSecondCam(:, 3), CurveSecondCam(:, 1:2), CurveFirstCam(:,3) , 'spline');

% triangulate
Trajectory = triangulate(CurveFirstCam(:, 1:2), MatchSecondCam, paramsFirstCam, paramsSecondCam);
%Trajectory = triangulate(CurveFirstCam(1:225, 1:2), CurveSecondCam(:, 1:2), paramsFirstCam, paramsSecondCam);

% plot trajectory
plotTrajectory(Trajectory);


function plotTrajectory(Trajectory)
figure;

curve = animatedline('lineWidth', 2);
set(gca, 'XLim', [-500 2800], 'YLim', [-300 3200], 'ZLim', [-400 1000]);
view(430, 240);
plotTable(-max(Trajectory(:,3)));
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


