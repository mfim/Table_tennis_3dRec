function cameraParams = cb_cab_r() 
    
    % Define images to process
    imageFileNames = {'./Calibration/R/6.png',...
        './Calibration/R/5.png',...
        './Calibration/R/3.png',...
        './Calibration/R/4.png',...
        './Calibration/R/2.png',...
        './Calibration/R/1.png',...
        };

        % Detect checkerboards in images
    [imagePoints, boardSize, imagesUsed] = detectCheckerboardPoints(imageFileNames);
    imageFileNames = imageFileNames(imagesUsed);

    % Read the first image to obtain image size
    originalImage = imread(imageFileNames{1});
    [mrows, ncols, ~] = size(originalImage);

    % Generate world coordinates of the corners of the squares
    squareSize = 23;  % in units of 'millimeters'
    worldPoints = generateCheckerboardPoints(boardSize, squareSize);

    % Calibrate the camera
    [cameraParams, imagesUsed, estimationErrors] = estimateCameraParameters(imagePoints, worldPoints, ...
        'EstimateSkew', false, 'EstimateTangentialDistortion', false, ...
        'NumRadialDistortionCoefficients', 2, 'WorldUnits', 'millimeters', ...
        'InitialIntrinsicMatrix', [], 'InitialRadialDistortion', [], ...
        'ImageSize', [mrows, ncols]);

    % Visualize pattern locations
    %h2=figure; showExtrinsics(cameraParams, 'CameraCentric');

    % Display parameter estimation errors
    displayErrors(estimationErrors, cameraParams);

end

