function [offset] = syncCam(BounceCoordFirstCam, BounceCoordSecondCam, BounceTsFirstCam, BounceTsSecondCam, paramsFirstCam, paramsSecondCam) 

% Calc real-world point of bounce
%  BounceCoordFirstCam(:,3) = 1;
%  BounceCoordSecondCam(:,3) = 1;
%  %RealWorldBounceFirstCam = pointsToWorld(paramsFirstCam, rotationMatrixFirstCam, translationVectorFirstCam, BounceCoordFirstCam);
%  %RealWorldBounceSecondCam = pointsToWorld(paramsSecondCam, rotationMatrixSecondCam, translationVectorSecondCam, BounceCoordSecondCam);
%  
%  RealWorldBounceFirstCam = BounceCoordFirstCam/paramsFirstCam;
%  RealWorldBounceFirstCam =  RealWorldBounceFirstCam./RealWorldBounceFirstCam(:,4); 
% 
%  RealWorldBounceSecondCam = BounceCoordSecondCam/paramsSecondCam;
%  RealWorldBounceSecondCam =  RealWorldBounceSecondCam./RealWorldBounceSecondCam(:,4); 
%BounceCoordFirstCam(:,4) = 1;
%BounceCoordSecondCam(:,4) = 1;
%RealWorldBounceFirstCam = BounceCoordFirstCam;
%RealWorldBounceSecondCam = BounceCoordSecondCam;

matches = inf([length(BounceCoordFirstCam(:,1)), 3]);


 for i = 1:length(BounceCoordFirstCam(:,1))
    for j = 1:length(BounceCoordSecondCam(:,1))
            %Check which is a better match in a time manner 
            new = triangulate(BounceCoordFirstCam(i, :), BounceCoordSecondCam(j,:), paramsFirstCam, paramsSecondCam);
            
            % How close the ball should be to the table 
            if(-new(:,3) < 300)
                % which is time difference is smaller
                if(abs(BounceTsSecondCam(j) - BounceTsFirstCam(i)) < abs(matches(i,1))) 
                    matches(i,:) = [BounceTsSecondCam(j) - BounceTsFirstCam(i), -new(:,3), j];
                end
           end
    end
 end
 
% Resolve the case where two points have the same match
% get the least percentage error to before 
%  for i = 1:length(BounceCoordSecondCam(:,1))
%      for j = 1:length(BounceCoordFirstCam(:,1))
%          if(matches(i,3) == matches(j,3))
%              if(abs
             
    
  % decide which camera is before 
 negative = 0;
 positive = 0;
  for i = 1:length(matches(:,1))
    if(matches(i, 1)<0)
        negative = negative + 1;
    else 
        positive = positive + 1;
    end
 end

if(positive > negative)
     offset = mean(matches(matches(:,1)>0, 1));
else
     offset = mean(matches(matches(:,1)<0, 1));
end
 
 %offset = matches(1,1);
end

 