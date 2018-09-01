function [offset] = syncCam(BounceCoordFirstCam, BounceCoordSecondCam, BounceTsFirstCam, BounceTsSecondCam, ...
    paramsFirstCam, rotationMatrixFirstCam, translationVectorFirstCam, ...
    paramsSecondCam, rotationMatrixSecondCam, translationVectorSecondCam) 

% Calc real-world point of bounce
  BounceCoordFirstCam(:,3) = 1;
  BounceCoordSecondCam(:,3) = 1;
%  %RealWorldBounceFirstCam = pointsToWorld(paramsFirstCam, rotationMatrixFirstCam, translationVectorFirstCam, BounceCoordFirstCam);
%  %RealWorldBounceSecondCam = pointsToWorld(paramsSecondCam, rotationMatrixSecondCam, translationVectorSecondCam, BounceCoordSecondCam);
%  
%  RealWorldBounceFirstCam = BounceCoordFirstCam/paramsFirstCam;
%  RealWorldBounceFirstCam =  RealWorldBounceFirstCam./RealWorldBounceFirstCam(:,4); 
% 
%  RealWorldBounceSecondCam = BounceCoordSecondCam/paramsSecondCam;
%  RealWorldBounceSecondCam =  RealWorldBounceSecondCam./RealWorldBounceSecondCam(:,4); 
BounceCoordFirstCam(:,4) = 1;
BounceCoordSecondCam(:,4) = 1;
RealWorldBounceFirstCam = BounceCoordFirstCam;
RealWorldBounceSecondCam = BounceCoordSecondCam;

matches = inf([length(RealWorldBounceFirstCam(:,1)), 2]);


 for i = 1:length(RealWorldBounceFirstCam(:,1))
    for j = 1:length(RealWorldBounceSecondCam(:,1))
       % Check if bounce is in the same plane 
       if(abs((RealWorldBounceSecondCam(j,4) - RealWorldBounceFirstCam(i,4))/RealWorldBounceFirstCam(i,4))< 0.2)
           %Check which is a better match in a time manner 
           if(abs(BounceTsSecondCam(j) - BounceTsFirstCam(i)) < abs(matches(i,1))) 
               matches(i,:) = [BounceTsSecondCam(j) - BounceTsFirstCam(i), j];
           end
       end
    end
 end
 
 offset = matches(1,1);
 
end

 