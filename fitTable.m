function [ H ] = fitTable( c1, c2, c3, c4 )
%Fit Ping Pong Table Homography 
%   Input: 4 image coordinates of the table ordered as:
%          left up - left down - right down - right up
%   Return the projective homography that bring the table 
%   to real world coordinates

    pi = [c1(1) , c1(2);  c2(1) , c2(2); c3(1) , c3(2); c4(1) , c4(2);];

    pt = [ 0,    0 ; 0,    1525; 2740, 1525 ;  2740, 0;];

     tform = fitgeotrans(pi, pt, 'projective');
     H = tform.T';
end

