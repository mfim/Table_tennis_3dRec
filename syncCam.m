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
            % can be fine-tuned 
            if(-new(:,3) < 300)
                % which time difference is smaller
                if(abs(BounceTsSecondCam(j) - BounceTsFirstCam(i)) < abs(matches(i,1))) 
                    matches(i,:) = [BounceTsSecondCam(j) - BounceTsFirstCam(i), -new(:,3), j];
                end
           end
    end
 end
 
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
     correct_matches = matches(matches(:,1)>0, :);
else
     correct_matches = matches(matches(:,1)<0, :);
end
 
% Resolve the case where two points have the same match
% get the least percentage error to before 
  for i = 1:length(correct_matches(:,1))
      for j = i:length(correct_matches(:,1))
          if((i ~= j) && correct_matches(i,3) == correct_matches(j,3))
              if(abs(correct_matches(i,3) - mean(correct_matches([1:i-1 i+1:j-1 j+1:end], 1))) ...
                      < abs(correct_matches(j,3) - mean(correct_matches([1:i-1 i+1:j-1 j+1:end], 1))))
                  final = correct_matches([1:i-1 i+1:end],:);
              else
                  final = correct_matches([1:j-1 j+1:end],:);
              end
          end
      end
  end
        
 offset = mean(final(:,1));
end

 