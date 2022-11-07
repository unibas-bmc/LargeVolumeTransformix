% read stack of tif slices and extract region of interest (ROI)
%
% [vol] = stackreader(fDir,roi,bgInt)
%
% fDir      tif filenames from fDir=dir([loadDir '*.tif']), done once to save time
% roi       image region, I(roi(3):roi(4),roi(1):roi(2),roi(5):roi(6))
% bgInt     background intensity, default minimum of datatype
%
% vol       3D array of extracted ROI 
%
function [vol] = stackreader(fDir,roi,bgInt)

%fDir = dir([loadDir '*.tif']);
if isempty(fDir)
    disp('Empty list of *.tif filenames')
    return
end
zlist = roi(5):roi(6);

info = imfinfo([fDir(1).folder filesep fDir(1).name]);

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
            'File: ' fDir(1).folder filesep fDir(1).name '\n' ...
            'BitsPerSample: ' num2str(info.BitsPerSample) '\n' ...
            'SampleFormat: ' info.SampleFormat '\n'];
    error(erstr)
end

% number of voxels
nv(1)=length(roi(3):roi(4));
nv(2)=length(roi(1):roi(2));
nv(3)=length(zlist);
% initialize with background intensity
if nargin<3
    bgInt=intmin(fmt);  
end
vol = bgInt*ones(nv(1),nv(2),nv(3),fmt);

% keep roi within image
newRoi(1)=max(1,roi(1));
newRoi(2)=min(info.Width,roi(2));
newRoi(3)=max(1,roi(3));
newRoi(4)=min(info.Height,roi(4));
newRoi(5)=max(1,roi(5));
newRoi(6)=min(length(fDir),roi(6));
newZlist=newRoi(5):newRoi(6);

diffRoi=roi-newRoi;
diffRoiAbs = sum(abs(diffRoi));  

if diffRoiAbs==0
    % roi is inside image
    for i = 1:nv(3)
        vol(:,:,i) = imread([fDir(zlist(i)).folder filesep fDir(zlist(i)).name],...
            'PixelRegion',{[roi(3), roi(4)],[roi(1), roi(2)]});
    end
else
    % roi is (partially) outside image
    if ~isempty(newZlist) && ~isempty(newRoi(3):newRoi(4)) && ~isempty(newRoi(1):newRoi(2))
        % none of the dimensions are empty
        for i = 1:length(newZlist)
            
            tmpI = imread([fDir(newZlist(i)).folder filesep fDir(newZlist(i)).name],...
                'PixelRegion',{[newRoi(3), newRoi(4)],[newRoi(1), newRoi(2)]});
            
            % put in right location
            vol(1-diffRoi(3):nv(1)-diffRoi(4),1-diffRoi(1):nv(2)-diffRoi(2),i-diffRoi(5))=tmpI;
        end
    end
end

end

