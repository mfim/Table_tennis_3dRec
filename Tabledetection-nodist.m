@clc;
clear;
close all;

%% Line Detection
I = imread('Images/tout.jpg');
WIDTH = max(size(I));

[lines, corners] = table_detection(I,1);
%
% [lines, corners] = table_detection(I,1);
lu = [586, 323, 1];
lm = [380, 371, 1];
ld = [67, 443, 1];

mu = [1054, 374, 1];
md = [645, 559, 1];

ru = [1729, 453, 1];
rm = [1702, 563, 1];
rd = [1656, 760, 1];

% Drawing
figure, imshow(I), hold on;

    text(lu(1),lu(2),'lu', 'Color', 'red');
    text(ld(1),ld(2),'ld', 'Color', 'red');
    text(ru(1),ru(2),'ru', 'Color', 'red');
    text(rd(1),rd(2),'rd', 'Color', 'red');

    text(md(1),md(2),'md', 'Color', 'red');
    text(mu(1),mu(2),'mu', 'Color', 'red');

    text(rm(1),rm(2),'rm', 'Color', 'red');
    text(lm(1),lm(2),'lm', 'Color', 'red');

drawline(ld, lu, 'red');
drawline(rd, ru, 'red');

drawline(ru, lu, 'yellow');
drawline(rd, ld, 'yellow');

drawline(mu, md, 'green');
drawline(lm, rm, 'green');

drawline(v1,v2,'red');

drawline(v1,lu,'red');
drawline(v1,ru,'red');

drawline(v2,lu,'red');
drawline(v2,ld,'red');

lline1 = cross(lu, ru);
lline2 = cross(ld, rd);
sline1 = cross(lu, ld);
sline2 = cross(ru, rd);

%lline1 = lline1/lline1(3);
%lline2 = lline2/lline2(3);
%sline1 = sline1/sline1(3);
%sline2 = sline2/sline2(3);

% Vanishing point
v1 = cross(sline1, sline2);
v2 = cross(lline1, lline2);

v1 = v1/v1(3);
v2 = v2/v2(3);
% Vanishing line
linf = cross(v1, v2);
linf = linf/linf(3);

% Affine transformation
H = [1, 0, 0; 0, 1, 0; linf];
Hs= [1, 0, 0; 0, 1, 0; linf];

T_aff = projective2d(Hs');

affine_im = imwarp(I, T_aff);
figure, imshow(affine_im);
title ('Affine reconstruction');
%%
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
