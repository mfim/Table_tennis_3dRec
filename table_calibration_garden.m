@clc;
clear;
close all;

I = imread('Images/j1.jpg');
WIDTH = max(size(I));

K = rgb2gray(I);
BW = edge(K,'Sobel');
[H,T,R] = hough(BW,'RhoResolution',0.5,'ThetaResolution',0.5);
P  = houghpeaks(H,30,'threshold',ceil(0.2*max(H(:))));
lines = houghlines(BW,T,R,P,'FillGap',300,'MinLength', 150);

figure, imshow(I), hold on;


for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',1,'Color', 'black');
   text(lines(k).point1(1), lines(k).point1(2), int2str(k), 'FontSize',15, 'Color', 'red'); 
end

lu = [1834, 473, 1];
lm = [1303, 669, 1];
ld = [469, 915, 1];

mu = [3012, 622, 1];
md =  [1851, 1290, 1];

ru = [4830, 827, 1];
rm = [4629, 1243, 1];
rd = [4391, 1949, 1];

text(lu(1),lu(2),'lu', 'Color', 'yellow');
text(ld(1),ld(2),'ld', 'Color', 'yellow');
text(ru(1),ru(2),'ru', 'Color', 'yellow');
text(rd(1),rd(2),'rd', 'Color', 'yellow');

text(md(1),md(2),'md', 'Color', 'yellow');
text(mu(1),mu(2),'mu', 'Color', 'yellow');

text(rm(1),rm(2),'rm', 'Color', 'yellow');
text(lm(1),lm(2),'lm', 'Color', 'yellow');

uline = points2line(lu, ru);
dline = points2line(ld, rd);
lline = points2line(lu, ld);
rline = points2line(ru, rd);
smline= points2line(mu, md);
slline= points2line(lm, rm);

table = struct('uline', points2line(lu, ru), ...  % upper line
               'dline', points2line(ld, rd), ...  % down line
               'lline', points2line(lu, ld), ...  % left line
               'rline', points2line(ru, rd), ...  % rignt line
               'smline', points2line(mu, md), ... % short middle line 
               'lmline', points2line(lm, rm) ); % long middle line
                
%               'lu', lu, 
%               'ld', ld, 'md', md, 'rd', rd, 'rm', rm, ... 
%               'ru', ru, 'lm', lm, 'mu'
           
lline1 = cross(lu, ru);
lline2 = cross(ld, rd);
sline1 = cross(lu, ld);
sline2 = cross(ru, rd);

v1 = cross(sline1, sline2);
v2 = cross(lline1, lline2);

v1 = v1/v1(3);
v2 = v2/v2(3);

% Vanishing line
linf = points2line(v1,v2);

%% Start of calibration

% pairs of orthogonal lines

ortho_pairs = [55, 60, 33, 42, 34, 39]; 

pairs = [];

figure, imshow(I), hold on;

for i = 1:2:length(ortho_pairs)
    
    il1 = ortho_pairs(i);
    il2 = ortho_pairs(i+1);
    
    l =  points2line(lines(il1).point1, lines(il1).point2);
    ol = points2line(lines(il2).point1, lines(il2).point2); 
   
    pairs = [ pairs, ...
              struct('l', l, 'ol', ol)]; 
    
    drawline(lines(il1).point1, lines(il1).point2, 'blue');
    drawline(lines(il2).point1, lines(il2).point2, 'blue');

end

% add table pairs 
pairs = [pairs, ... 
         struct('l', table.uline, 'ol', table.lline), ...  
         struct('l', table.dline, 'ol', table.rline), ... 
         struct('l', table.lline, 'ol', table.lmline), ... 
         struct('l', table.dline, 'ol', table.smline) ];

         %struct('l', table.uline, 'ol', table.rline), ...  
         %struct('l', table.dline, 'ol', table.lline), ...  
     
% Vanishing point pairs 

vppairs = [];
     
for i=1:length(pairs) 
    
    vp = points2line(pairs(i).l,  linf);          % Duality 
    vpo = points2line(pairs(i).ol, linf);          %
    
    vppairs = [ vppairs, ... 
                struct('vp', vp, 'vpo', vpo )];
end 

% Check vppairs it is inconstitently constrained 
% A lot of vanishing points are doubled

