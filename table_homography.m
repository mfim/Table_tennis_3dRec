@clc;
clear;
close all;

I = imread('ping.png');
WIDTH = max(size(I));

[lines, corners] = table_detection(I, 1);

% ping;

for k = 1:length(lines)
 
  a = lines(k).point1; 
  b = lines(k).point2;
  plot([a(1), b(1)], [a(2), b(2)],'LineWidth',2,'Color', 'black');
  text(b(1), b(2), num2str(k), 'Color', 'red');  
  
end 

% Middle Lines

i_mll = 2; 
i_msl = 6; 

mlline = points2line(lines(i_mll).point1, lines(i_mll).point2);
msline = points2line(lines(i_msl).point1, lines(i_msl).point2);

ru = find(corners(:,1) == max(corners(:,1)));
ld = find(corners(:,1) == min(corners(:,1)));
lu = find(corners(:,2) == min(corners(:,2)));
rd = find(corners(:,2) == max(corners(:,2)));
% Homogeneous coordinates for corners
c1 = corners(lu,:); c2 = corners(ld,:); c3 = corners(rd,:); c4 = corners(ru,:);

H = fitTable(c1, c2, c3, c4);
  
Tf = projective2d(H');
reconstructed_im = imwarp(I, Tf);
figure, imshow(reconstructed_im);
title('Shape reconstruction');
 
 
 