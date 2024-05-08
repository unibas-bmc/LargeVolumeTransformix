% Binning for image registration
%
% does in parallel per binned slice
% - loading input tif slices needed for binning output
% - changes intensity to the specified range
% - writes binned tif slice
%
% Note:
%
% afterwards to save as .mha file:
% import binned tiff stack in FIJI/ImageJ, 
% z-direction is flipped (image->transform->flip Z)
% voxel size set (image->properties) and saved as mha
%
% if origin of b32 is [0 0 0], 
% then origin of b16 should be [-res -res -res], 
% origin can be set in FIJI, BUT is not take over when saving mha!!!
%
bf1 = 8;   % bin factor of input tiff images

bf2 = 32;  % bin factor of output tiff images

this_bf = floor(bf2/bf1);

% data directory, e.g. link this to absolute path
dataDir = ['dataDir' filesep];

if bf1>1
    loadDir = [dataDir 'reco_b' num2str(bf1) filesep];
else
    loadDir = [dataDir 'reco' filesep];
end
thisDir = [dataDir 'reco_b' num2str(bf2) '_forRegi' filesep];

if not(isfolder(thisDir)); mkdir(thisDir); end

dirInfo = dir([loadDir 'reco_*.tif']);
numIms = length(dirInfo);
im = imread([dirInfo(1).folder filesep dirInfo(1).name]);
[sy,sx] = size(im);
if bf1 == 0
    osy = this_bf*floor((sy-1)/this_bf); osx = this_bf*floor((sx-1)/this_bf);
else
    osy = this_bf*floor(sy/this_bf); osx = this_bf*floor(sx/this_bf);
end

numImsFinal = floor(numIms/this_bf);

tagstruc.ImageLength = osy/this_bf; % y
tagstruc.ImageWidth = osx/this_bf; % x
tagstruc.BitsPerSample = 16;
tagstruc.SamplesPerPixel = 1;
tagstruc.Compression = Tiff.Compression.None;
tagstruc.SampleFormat = Tiff.SampleFormat.UInt;
tagstruc.Photometric = Tiff.Photometric.MinIsBlack;
tagstruc.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

% intensity scaling
gs = [-20000,-10000];  % max and min intensity
disp(['Image intensities are scaled from range [' num2str(gs(1)) ',' num2str(gs(2)) '] to uint16'])

parfor b = 1:numImsFinal
    outputFile = [thisDir 'reco_b' num2str(bf2) '_' num2str(b,'%04d') '.tif'];
    if exist(outputFile,'file')
        disp([outputFile ' exists'])
    else
        t00 = tic;
        
        i1 = ((b-1)*this_bf) + 1;
        
        vol = zeros(osy,osx,this_bf,'int16');
        for i = 1:this_bf
            im = imread([dirInfo(i1-1+i).folder filesep dirInfo(i1-1+i).name]);
            vol(:,:,i) = im(1:osy,1:osx);
        end
        volb = reshape(double(vol),this_bf,osy/this_bf,this_bf,osx/this_bf,this_bf,1);
        volb = squeeze(mean(volb,[1,3,5]));        % take mean over this_bf^3 voxels
        volb = (2^16)*(volb-gs(1))/(gs(2)-gs(1));  % intensity change
        
        t = Tiff(outputFile, 'w');
        t.setTag(tagstruc); t.write(uint16(volb)); t.close();
        
        fprintf('Binned block %i of %i \n',b,numImsFinal)
        toc(t00)
    end
end


