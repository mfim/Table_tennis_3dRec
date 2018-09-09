function line = points2line(p1, p2)
% Function that takes the image coordinates of 2 points and returns 
% the homogeneous coordinates of the line that goes through these two
% points
  
  p1 = [p1(1), p1(2), 1];
  p2 = [p2(1), p2(2), 1];
  
  line = cross(p1, p2); 
  line = line/line(3);

end

