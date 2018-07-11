@clc;
clear;
close all;

%% Edging
I = imread('ping2.png');
[HEIGHT, WIDTH, DIM] = size(I);

LAB = rgb2lab(I);
HSV = rgb2hsv(I);

rLAB = rangefilt(LAB);
 
rLAB = rLAB(:,:,1);

K = rgb2gray(I);
BW = edge(K,'Sobel');
% figure; imshow(BW); 

%% Line Detection 

[H,T,R] = hough(BW,'RhoResolution',0.5,'ThetaResolution',0.5);

P  = houghpeaks(H,20,'threshold',ceil(0.2*max(H(:))));
x = T(P(:,2)); 
y = R(P(:,1));

plotInter = 0; 

% Find lines and plot them
lines = houghlines(BW,T,R,P,'FillGap',150,'MinLength', 300);
figure, imshow(I), hold on;

fLines = [];

couleurs = ['r','b','y','g']; 

COLOR_TRESHOLD = 40;
BLUE_MAX = 270; 
BLUE_MIN = 200;
VLIGHT_MIN = 20;
MIN_AREA = 1e+4; 
MAX_AREA = 2*1e+4; 

for k = 1:length(lines)
 
  a = lines(k).point1; 
  b = lines(k).point2;
  
  coefficients = polyfit([a(1), b(1)], [a(2), b(2)], 1);
  
  lines(k).a = coefficients(1); 
  lines(k).b = coefficients(2);
  
  lines(k).centroid = (lines(k).point1 + lines(k).point2)/2;
  
  centroid = lines(k).centroid;
  
  % Orthogonal line that passes through the center 
  orthoa = -1/coefficients(1); 
  orthob = -orthoa*centroid(1) + centroid(2) ;
  % Orthogonal line equation  that passes through the center 
  ortho = @(x) orthoa*x + orthob;
  
  %Plot orthogonal linesf
  % x=[0,width];
  % plot(x, ortho(x), 'Color','cyan', 'LineWidth',2);   
  
  shift = 20 / sqrt(1+orthoa^2); 
  % +1 because the HSC doesn't like zeros
  s1 = [ceil(centroid(1) + shift), ceil(ortho(centroid(1) + shift))+1];
  s2 = [ceil(centroid(1) - shift), ceil(ortho(centroid(1) - shift))+1];
  
  if plotInter
    plot(s1(1), s1(2) ,'o','LineWidth',2,'Color','yellow');
    plot(s2(1), s2(2) ,'o','LineWidth',2,'Color','yellow');
  end

  blue1 = HSV(s1(2), s1(1), 1)*360;
  blue1V = HSV(s1(2), s1(1), 3)*100;
  
  blue2 = HSV(s2(2), s2(1), 1)*360;
  blue2V = HSV(s2(2), s2(1), 3)*100;
  
  lines(k).blue1 = blue1 > BLUE_MIN & blue1 < BLUE_MAX & blue1V > VLIGHT_MIN;
  lines(k).blue2 = blue2 > BLUE_MIN & blue2 < BLUE_MAX & blue2V > VLIGHT_MIN;
  
  lines(k).bluediff = xor(lines(k).blue1, lines(k).blue2);
  
  % Filter on the line orientation and if it is a separating line 
  if abs(lines(k).theta) > 15 & lines(k).bluediff
      fLines = [fLines; lines(k)];
  end      
end

allLines = lines; 

lines = fLines;

% Hierarchial clustering on the table border lines

A = [lines.a ; lines.b]';

 %distances = squareform(pdist(A));

 % linkages  = linkage(A, 'complete', 'seuclidean');
 % linesCluster = cluster(linkages, 'maxclust', 4);


area = 0; 
iterations = 0; 
areaArray = []
while (iterations < 1)
    
    linesCluster = kmeans(A,4,'Distance','sqeuclidean');
    for k = 1:length(fLines)
  
    xy = [fLines(k).point1; fLines(k).point2];
    disp(xy);
  
    % Assign cluster in line list
    fLines(k).cluster = linesCluster(k);
  
    fLines(k).length = sqrt( (xy(2,1)-xy(1,1))^2 + (xy(2,2)-x(1,2))^2 );
  
     if plotInter
        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color', couleurs(linesCluster(k)));
        %Show Theta of the lines 
        text(xy(1,1),xy(1,2),num2str(lines(k).theta), 'Color', 'red');
        % plot beginnings and ends of lines
        plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
        plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');   
     end
    end

    % Find horizontal lines of the table based on approximative length

    table = [];

    for c=1:max(linesCluster)
    
        % Lines are defined with ax + b 
        % so we are trying to find those for our new line 
        % which is a weighted average of all of these
    
        currentClust = fLines([fLines.cluster] == c);
        % This is to make longer vectors weigh more
        vlength = [currentClust.length].^4; 
        % Calculate total length for norm 
        total = sum(vlength);
        % Divide by the total to get the weights of each line
        weights  = vlength/total; 
    
        % Calculate weighted coefficients 
        coeff.a = dot([currentClust.a]', weights);
        coeff.b = dot([currentClust.b]', weights);
    
        table = [table; coeff];  
    end
    
    icpt = [];
    icptpoints = []; 
    
    % Calculate intercepts
    for i=1:length(table)
  
        f1 = @(x) x*table(i).a + table(i).b; 
   
        x =  [0, WIDTH]; % linspace(0, width); % Need only two point to plot
    
        for j = i+1:length(table)
            f2 = @(x) x*table(j).a + table(j).b; 
            icpt(i, j) = fzero(@(x) f1(x) - f2(x), 1*((-1)^j));      
            icptpoints = [icptpoints ; [icpt(i, j), f1(icpt(i, j))]];
            % plot(icpt(i, j), f1(icpt(i, j)),'o', 'Color','red', 'LineWidth', 2); 
            
        end      
    end
    
    corners = icptpoints(icptpoints(:,1) > 0 & icptpoints(:,2) > 0, :);
    area = polyarea(corners(:,1), corners(:,2));
    iterations = iterations + 1;
    areaArray = [areaArray; area]
end 
% End of loop 

%Plotting final lines
for i=1:length(table)
    
    plot(x, f1(x), 'Color','blue', 'LineWidth',2); 
    f1 = @(x) x*table(i).a + table(i).b; 
   
    x =  [0, WIDTH];% Need only two point to plot
        
    for j = i+1:length(table)
        f2 = @(x) x*table(j).a + table(j).b; 
        temp = fzero(@(x) f1(x) - f2(x), 1*((-1)^j));      
        plot(temp, f1(temp),'o', 'Color','red', 'LineWidth', 2);     
    end  
end 

%    







            


