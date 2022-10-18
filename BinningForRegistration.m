%% Binning for registration
bf1 = 8;
bf2 = 32;

this_bf = floor(bf2/bf1);

loadDir = ['/media/griffin/anatomix_mousebrain/mouse4_perf_eth/reco_b' num2str(bf1) filesep];
saveDir = '/media/griffin/anatomix_mousebrain/mouse4_perf_eth/';
thisDir = [saveDir 'reco_b' num2str(bf2) '_forRegi' filesep];
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

gs = [-20000,-10000];

parfor b = 1:numImsFinal
    t00 = tic;
    
    i1 = ((b-1)*this_bf) + 1;
    
    vol = zeros(osy,osx,this_bf,'int16');
    for i = 1:this_bf
        im = imread([dirInfo(i1-1+i).folder filesep dirInfo(i1-1+i).name]);
        vol(:,:,i) = im(1:osy,1:osx);
    end
    volb = reshape(double(vol),this_bf,osy/this_bf,this_bf,osx/this_bf,this_bf,1);
    volb = squeeze(mean(volb,[1,3,5]));
    volb = (2^16)*(volb-gs(1))/(gs(2)-gs(1));
    
    t = Tiff([thisDir 'reco_b' num2str(bf2) '_' num2str(b,'%04d') '.tif'], 'w');
    t.setTag(tagstruc); t.write(uint16(volb)); t.close();
    
    fprintf('Binned block %i of %i \n',b,numImsFinal)
    toc(t00)
end








