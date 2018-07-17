clc;
clear;
close all;
fontSize = 20;

MAX_ITERATIONS = 3;
BALL_SIZE = 22; % pixels of the ball 
% just for now.. keep the ballcolor
ballColor = 'o';
first_threshold = 3; second_threshold = 8;

firstPosition = 'Video/5-m-60fps.mp4';
%firstPosition = 'Video/5-m-120fps.mp4';
secondPosition = 'Video/5-r-60fps.mp4';
%secondPosition = 'Video/5-r-120fps.mp4';
%videoName = 'Video/1-m-30fps.mp4';


% prompt = 'Frames to be skipped: ';
% skips = str2double(input(prompt, 's'));
% if isnan(skips) || fix(skips) ~= skips
%   disp('Please enter an integer');
%   return;
% end

[PositionsFirstCam, CurveFirstCam, BounceTsFirstCam, BounceCoordFirstCam, StrikeTsFirstCam, StrikeCoordFirstCam] = ballTracking(firstPosition, ... 
     MAX_ITERATIONS, BALL_SIZE, ballColor, first_threshold, second_threshold, 30);
[PositionsSecondCam, CurveSecondCam, BounceTsSecondCam, BounceCoordSecondCam, StrikeTsSecondCam, StrikeCoordSecondCam] = ballTracking(secondPosition,...
     MAX_ITERATIONS, BALL_SIZE, ballColor, first_threshold, second_threshold, 0);

 % Bounce 
 % if only one that's the right
 % if more compare time diff to match
 % sep bounce from stroke 
 
% calc delta from two videos 
% set the times right 
% plot trajectory from two 
 
 
% [xData, yData, zData] = prepareSurfaceData(time, x, y );
% 
% % Set up fittype and options.
% ft = 'cubicinterp';
% 
% % Fit model to data.
% [fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );
% 
% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, [xData, yData], zData );
% legend( h, 'untitled fit 1', 'time vs. x, y', 'Location', 'NorthEast' );
% % Label axes
% xlabel x
% ylabel y
% zlabel time
% grid on
% view( 0.0, 90.0 );
%  
