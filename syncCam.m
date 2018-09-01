function [offset] = syncCam(BounceCoordFirstCam, BounceCoordSecondCam, BounceTsFirstCam, BounceTsSecondCam) 

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

matches = inf([length(BounceCoordFirstCam(:,1)), 2]);


 for i = 1:length(BounceCoordFirstCam(:,1))
    for j = 1:length(BounceCoordSecondCam(:,1))
            %Check which is a better match in a time manner 
           if(abs(BounceTsSecondCam(j) - BounceTsFirstCam(i)) < abs(matches(i,1))) 
               matches(i,:) = [BounceTsSecondCam(j) - BounceTsFirstCam(i), j];
           end
    end
 end
 
 
 %offset = matches(1,1);
 offset = mean(matches(:,1));
end

 