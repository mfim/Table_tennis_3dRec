I = imread('ping2.png');

c = [1168.26748730376;822.692259078109;343.946154724962;930.826492215078;];
r = [132.307612787839;77.6503062827622;381.782977945311;567.312891087348];

BW = roipoly(I,c,r);

figure
imshow(BW)