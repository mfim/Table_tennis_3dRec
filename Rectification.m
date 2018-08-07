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
c1=[c1,1]; c2=[c2,1]; c3 =[c3,1]; c4=[c4,1]; 
% Plot corners 
text(c1(1),c1(2),'C1', 'Color', 'yellow');
text(c2(1),c2(2),'C2', 'Color', 'yellow');
text(c3(1),c3(2),'C3', 'Color', 'yellow');
text(c4(1),c4(2),'C4', 'Color', 'yellow');
% Find edge lines 
sline1 = points2line(c1, c2); 
sline2 = points2line(c3, c4); 
lline1 = points2line(c1, c4); 
lline2 = points2line(c2, c3);

% Vanishing point
v1 = cross(sline1, sline2); 
v2 = cross(lline1, lline2);
% 
v1 = v1/v1(3);
v2 = v2/v2(3);
% Vanishing line
linf = cross(v1, v2);
linf = linf/linf(3);
% Affine transformation
Haff = [1, 0, 0; 0, 1, 0; linf];
T_aff = projective2d(Haff');
affine_im = imwarp(I, T_aff);
figure, imshow(affine_im);
title ('Affine reconstruction');

% Each line of this matrice contains two orthogonal linesx
ortholines = [  sline2/Haff, lline2/Haff; ... 
                sline1/Haff, lline1/Haff; ...
                sline2/Haff , mlline/Haff;  ... 
                lline1/Haff , msline/Haff; 
              ];

% Least squares to find dual degenerate conic 
X = zeros(4,2);
t = zeros(4,1);

for i=1:4
   
    X(i,1) = ortholines(i,1)*ortholines(i,4);
    X(i,2) = ortholines(i,1)*ortholines(i,5) + ortholines(i,2)*ortholines(i,4);
    
    t(i) = - ortholines(i,2) * ortholines(i,5);
       
end

w = X\t;

C = [w(1), w(2), 0;
     w(2),    1, 0;
     0,       0, 0];
 
[U, S, V] = svd(C);

Ht = U * diag([sqrt(S(1, 1)), sqrt(S(2, 2)), 1]);

Hrecon = Ht \ Haff;

T = projective2d(Hrecon');
reconstructed_im = imwarp(I, T);
figure, imshow(reconstructed_im);
title('Shape reconstruction');

cc1 = c1 / Hrecon;
cc2 = c2 / Hrecon; 
cc4 = c4 / Hrecon; 

cc1 = cc1 / cc1(3); 
cc2 = cc2 / cc2(3);
cc4 = cc4 / cc4(3);

long = norm(cc1(1:2) - cc4(1:2), 2); 
short= norm(cc1(1:2) - cc2(1:2), 2); 

ratio = long / short; 

if ratio > 1
    ratio = 1/ratio
end





