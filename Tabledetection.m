@clc;
clear;
close all;

%% Line Detection 
I = imread('ping2.png');
WIDTH = max(size(I));
% 
[lines, corners] = table_detection(I,1);

%% AFFINE RECONSTRUCTION
% Line Normalisation

% Reorganisation of the interesect points (corners)
ru = find(corners(:,1) == max(corners(:,1)));
ld = find(corners(:,1) == min(corners(:,1)));
lu = find(corners(:,2) == min(corners(:,2)));
rd = find(corners(:,2) == max(corners(:,2)));

c1 = corners(lu,:); c2 = corners(ld,:); c3 = corners(rd,:); c4 = corners(ru,:);
c1=[c1,1]; c2=[c2,1]; c3 =[c3,1]; c4=[c4,1]; 

text(c1(1),c1(2),'C1', 'Color', 'yellow');
text(c2(1),c2(2),'C2', 'Color', 'yellow');
text(c3(1),c3(2),'C3', 'Color', 'yellow');
text(c4(1),c4(2),'C4', 'Color', 'yellow');

lline1 = cross(c1, c4); 
lline2 = cross(c2, c3); 
sline1 = cross(c1, c2); 
sline2 = cross(c3, c4);

lline1 = lline1/lline1(3); 
lline2 = lline2/lline2(3); 
sline1 = sline1/sline1(3); 
sline2 = sline2/sline2(3); 

N = diag([WIDTH, WIDTH, 1]);

nlline1 = lline1/N; 
nlline2 = lline2/N; 
nsline1 = sline1/N; 
nsline2 = sline2/N;

% Vanishing point
v1 = cross(sline1/N, sline2/N); 
v2 = cross(lline1/N, lline2/N);
% Vanishing line
linf = cross(v1, v2);
linf = linf/linf(3);

% Affine transformation
H = [1, 0, 0; 0, 1, 0; linf];
Hs= [1, 0, 0; 0, 1, 0; linf*N];

T_aff = projective2d(Hs');

affine_im = imwarp(I, T_aff);
figure, imshow(affine_im);
title ('Affine reconstruction');

% Each line of this matrice contains two orthogonal linesx
ortholines = [nsline1/H, nlline2/H; ... 
              nsline2/H, nlline1/H; ... 
              nsline2/H, nlline2/H; ... 
              nsline1/H, nlline1/H];

% Least squares to find dual degenerate conic 
X = zeros(2,2);
t = zeros(2,1);

for i=1:2
   
    X(i,1) = ortholines(i,1)*ortholines(i,4);
    X(i,2) = ortholines(i,1)*ortholines(i,5) + ortholines(i,2)*ortholines(i,4);
    
    t(i) = - ortholines(i,2) * ortholines(i,5);
       
end

w = X\t;

C = [w(1), w(2), 0;
     w(2),    1, 0;
     0,       0, 0];
 
% https://stackoverflow.com/questions/12449695/metric-rectification-using-the-dual-degenerate-conic-in-matlab
 
[U, S, V] = svd(C);

sqrt = U * diag([sqrt(S(1, 1)), sqrt(S(2, 2)), 1]);

Hrecon = sqrt \ Hs;

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






