function [] = drawline(p1,p2,color)
%DRAWLINE Summary of this function goes here
%   Detailed explanation goes here

    plot([p1(1), p2(1)],[p1(2), p2(2)] ,'LineWidth',2,'Color', color);
    
end

