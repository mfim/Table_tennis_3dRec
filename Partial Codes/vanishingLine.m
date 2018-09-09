function [ linf ] = vanishingLine(lu, ld, rd, ru)
%VANISHINGLINE Given a table corners, this function returns the vanishing
%lines
%   Detailed explanation goes here
    lline1 = cross(lu, ru);
    lline2 = cross(ld, rd);
    sline1 = cross(lu, ld);
    sline2 = cross(ru, rd);

    % Vanishing point
    v1 = cross(sline1, sline2);
    v2 = cross(lline1, lline2);

    v1 = v1/v1(3);
    v2 = v2/v2(3);
    % Vanishing line
    linf = cross(v1, v2);
    linf = linf/linf(3);

end

