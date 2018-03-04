% Programming By ALireza Fasih
% Email: ar_fasih@yahoo.com
% Please copy images (*.ppm) to work directory, then push F5 to view results.

clc;       % clear consul
clear all; % clear all variable

% Ball Threshold Min & Max Values
th_min_ball=80;
th_max_ball=214;

% number of frams
frams = 88;

% Ball Pos initialize
old_ball_pos_x=0;
old_ball_pos_y=0;
Ball_speed_list(2,10)=0;
ball_vel_x=0;
ball_vel_y=0;

delta_x=25;
delta_y=25;


counter=0;
for i = 21: frams
 fname=strcat('images.',int2str(i),'.ppm'); % Incrimental String Variable for Indexing file name

 
 if (i<35)
  delta_x=35;
  delta_y=35;
 else
  delta_x=25;
  delta_y=25;     
 end
 
 if (i>80)
  delta_x=30;
  delta_y=15;
 end
 
 
 k=imread(fname);       % Load Image to 'k' Matrix 
 k1=k;
 if (old_ball_pos_x ~=0)
   k1=imcrop(k,[ round(old_ball_pos_x-delta_x) round(old_ball_pos_y-delta_y) delta_x.*2 delta_y.*2 ]);
 end
 
 I=rgb2gray(k1);        % Converting RGB Image to Gray Scale Image
 I=im2double(I);        % Converting Gray scale Image to Double type

  J = medfilt2(I,[3 3]); % Median Filter , 3x3 Convolution on Image 
 I2 = imadjust(J);       % Improve to quality of Image and adjusting contrast and brightness values
 Ib = I2> 0.8627;        %(220./255); % Binary Threshold for Ball
 

 [labeled,numObjects] = bwlabel(Ib,4);              % Indexing segments by binary label function
 graindata = regionprops(labeled,'all');            % extracting properties from regions, 'all' means we need to all attributes
 display(i); % show cycle

  %---------------------- Ball Area,  Filter --------------- 
  a_f = find( [graindata.Area] > th_min_ball & [graindata.Area] < th_max_ball); % Because Ball area is in range of 80 and 211  
  c1=length(a_f);
  %---------------------- Ball W/H Ratio,  Filter --------------- 
  sel=1;     % default value
  for j = 1 : c1   
    
    Ma=graindata(a_f(j)).MajorAxisLength;
    Mi=graindata(a_f(j)).MinorAxisLength;
    ratio_w_h =   Ma./Mi;
    
     if (ratio_w_h<1.7)  % Ball or semi circle ratio range = ratio>0.6 & ratio<1.7  
        sel = a_f(j);
     end
  end
  

 imshow(k); %k
 hold on;
 
 set(findobj(gca,'Type','line','Color',[0 0 1]),...
    'Color','Blue',...
    'LineWidth',3)

 % Draw point on Ball Center
 c1=graindata(sel).Centroid;

 pos_list(1,i+1)=c1(1);
 pos_list(2,i+1)=240-c1(2);
 plot(c1(1)+old_ball_pos_x-delta_x,c1(2)+old_ball_pos_y-delta_y,'rx');
 
 for tt=0:0.05:(2*pi)
  plot(c1(1)+old_ball_pos_x-delta_x+40*cos(tt),c1(2)+old_ball_pos_y-delta_y+40*sin(tt),'g');
 end
 

  if (i>0)
    ball_vel_x = (old_ball_pos_x - c1(1))./24; % 24 frame in second    velocity = variation of position in delta time    
    ball_vel_y = (old_ball_pos_y - c1(2))./24;
  end
    
  Ball_speed_list(1,i+1)=ball_vel_x;
  Ball_speed_list(2,i+1)=ball_vel_y;


 if (counter>0)
  old_ball_pos_x=old_ball_pos_x-delta_x+c1(1);
  old_ball_pos_y=old_ball_pos_y-delta_y+c1(2);  
 else
  old_ball_pos_x=c1(1);
  old_ball_pos_y=c1(2);       
 end
 counter=counter+1;
 
 pause(0.1);
 
end            

