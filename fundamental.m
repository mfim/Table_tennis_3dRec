firstPosition = 'Video/5-m-120fps.mp4';
secondPosition = 'Video/5-r-120fps.mp4';


v1 = VideoReader(firstPosition, 'CurrentTime', 1.5);
v2 = VideoReader(secondPosition, 'CurrentTime', 1.5);


view1 = readFrame(v1);
view2 = readFrame(v2);

% imshow(view1);

matchedPoints1= [358,379;
                 433,218;
                 923,571;
                 1035,500;
                 1104,229;
                 1067,364; 
                 1163,134;
                 819,76];

matchedPoints2= [447,154;
                 787,174;
                 121,280;
                 441,606;
                 397,330;
                 475,512;
                 971,572;
                 1151,292];

% figure; showMatchedFeatures(view1,view2,matchedPoints1,matchedPoints2);

F = estimateFundamentalMatrix(matchedPoints1,matchedPoints2, 'Method', 'Norm8Point');


x1 = [matchedPoints1(3,:), 1];
x2 = [matchedPoints2(3,:), 1]';