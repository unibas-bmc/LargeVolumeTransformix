function [sx,sy,sz] = GetStackSize(dirinfo)
%loaddir = [readDir '*.tif'];
%dirinfo = dir(loaddir);

%tmpim = imread([dirinfo(1).folder filesep dirinfo(1).name]);
%[sy,sx] = size(tmpim);
t = Tiff([dirinfo(1).folder filesep dirinfo(1).name]);
sx = getTag(t,'ImageWidth');
sy = getTag(t,'ImageLength');

sz = length(dirinfo);
end
