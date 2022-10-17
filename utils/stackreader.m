function [vol] = stackreader(loadDir,roi)

tmp = dir([loadDir '*.tif']);
zlist = roi(5):roi(6);

info = imfinfo([tmp(1).folder filesep tmp(1).name]);

if info.BitsPerSample == 32
    fmt = 'single';
elseif info.BitsPerSample == 64
    fmt = 'double';
elseif info.BitsPerSample == 16
    if strcmp(info.SampleFormat,"Two's complement signed integer")
        fmt = 'int16';
    else
        fmt = 'uint16';
    end
else
    erstr = ['Could not identify the data type of the image stack: \n' ...
            'File: ' tmp(1).folder filesep tmp(1).name '\n' ...
            'BitsPerSample: ' num2str(info.BitsPerSample) '\n' ...
            'SampleFormat: ' info.SampleFormat '\n'];
    error(erstr)
end

vol = zeros(length(roi(3):roi(4)),length(roi(1):roi(2)),length(zlist),fmt);
for i = 1:length(zlist)
    vol(:,:,i) = imread([tmp(zlist(i)).folder filesep tmp(zlist(i)).name],...
        'PixelRegion',{[roi(3), roi(4)],[roi(1), roi(2)]});
end

end

