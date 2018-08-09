
I1 = rgb2gray(imread('Images/j4.jpg'));
I2 = rgb2gray(imread('Images/j5.jpg'));

trsh = 1200;

imshow([I1,I2]);

points1 = detectSURFFeatures(I1, 'MetricThreshold', trsh);
points2 = detectSURFFeatures(I2, 'MetricThreshold', trsh);

[features1,valid_points1] = extractFeatures(I1,points1);
[features2,valid_points2] = extractFeatures(I2,points2);

indexPairs = matchFeatures(features1,features2);

matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);

figure; showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2);


fRANSAC = estimateFundamentalMatrix(matchedPoints1,...
    matchedPoints2,'Method','RANSAC',...
    'NumTrials',2000,'DistanceThreshold',1e-4)
