function [roi] = CoordsToBoundingboxRoi(coords)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

tmp = min(coords,[],1);
x1 = floor(tmp(1)); y1 = floor(tmp(2)); z1 = floor(tmp(3));
tmp = max(coords,[],1);
x2 = ceil(tmp(1)); y2 = ceil(tmp(2)); z2 = ceil(tmp(3));

roi = [x1,x2,y1,y2,z1,z2];

end

