function [cameraParameters, rotationMatrix, translationVector] = calcExtrinsics(frame, intrinsics)

% official table dimension in cm
worldPoints = [0, 0; 1525, 0; 0, 2740; 1525, 2740];

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

figure; plotCamera('Location', location, 'Orientation', orientation, 'Size', 200);
hold on;
grid on;    
pcshow([worldPoints, zeros(size(worldPoints, 1),1)], ...
    'VerticalAxisDir','down','MarkerSize',40);

end

