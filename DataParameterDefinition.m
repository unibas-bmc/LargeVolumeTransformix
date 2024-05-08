% Definition of parameters and data directories for TransformToAtlas.m
%
%
doTesting = 0;    % 1: check stitched result against warped full image
showROIBox = 0;   % 1: show bounding box in unwarped space during testing
isotropicSplit=0; % 1: same number of grids in each direction, else force 1x1xZ 
optimizedSplit=0; % 1: optimize split to have least slices, else force 1x1xZ

%% Pixel sizes and resolutions used for registration
ps_fac = 100; % divide all values by this when using itk coordinates
% this helps with registration, otherwise bending energy penalty is below
% precision and won't work

ps_atlas_10 = 10;
ps_atlas_25 = 25;
ps_uct_1x = 0.65; % [um]

bf_regi = 32; % bin factor for uct used in registration

ps_atlas_regi = ps_atlas_25;
ps_uct_regi = ps_uct_1x*bf_regi;

% which HR resolutions you want (i.e. binning factor)
% e.g. 32, 16, 8, 4, 2, 1
% do one resolution for all regions
res = 16;  

% increase roi of unwarped image to avoid interpolation artifacts
% define number of pixels to extend on each side, i.e. 5 pixels in bin32
extendBy = 5*round(32/res);

onlyFg = 1;  % transform only rois containing forground

% break into rois no larger than MaxSizeGB
if doTesting
    % TEST with smaller max size to force breaking up image
    MaxSizeGB = 0.2;  % splits 32x into 2x2x2
    testingStr = '_test';
else
    if res==2 || res==1
        isotropicSplit=1;  % unwarped skew, so bad with slices
    end
    if isotropicSplit
        MaxSizeGB = 6;  % 11 GB is swapping after half of the slices on my machine
    else
        % 10.1 GB unwarped and 8.5 GB warped was still OK (8x)
        MaxSizeGB = 11; 
    end
    testingStr = '';
end

%% Directory structure
% base directory where high resolution tiff stacks are stored
if res>4
    loadBaseDir = '/media/christine/soleil_july21_mouse4_eth/LargeVolumeTransformixData/reconstructions/mouse4_eth/';
else
    loadBaseDir = '/media/christine/anatomix_mousebrain/mouse4_perf_eth/';
end
% for bin 1, tiffs should be here: [loadBaseDir 'reco' filesep];
% for any other binning, tiffs should be here: [loadBaseDir 'reco_bX' filesep];
warpManual = 0;
if warpManual
    % this is the manual affine
    regStr = 'ManualAffine';
    % location of registration files
    regiBaseDir = ['.' filesep 'example' filesep 'registration' filesep];
    
    % base directory where warped images will be stored
    warpBaseDir = [regiBaseDir 'warping' regStr testingStr filesep];
    if not(isfolder(warpBaseDir)); mkdir(warpBaseDir); end
    
    % TransformParameters of low-res registration to apply to high res data
    regiDir = [regiBaseDir 'affine0/']
    lrRegistrationTransform = [regiDir 'TransformParameters.0.txt'];
    % low-resolution wared image
    lr_warped_filename = [regiBaseDir 'volumes/result.mha'];
else
    % this is the best non-rigid registration result
    % /media/christine/soleil_july21_mouse4_eth/atlas_registration
    % v6/non-rigid/testsGrad/Def_RegSeries_mask_LMall_gs16_bep1000_wLM0p1/result.0.mha 
    regStr = 'NR';
     
    % location of registration files
    regiBaseDir = '/media/christine/soleil_july21_mouse4_eth/atlas_registration/';
    
    % base directory where warped images will be stored
    warpBaseDir = [regiBaseDir 'warping' regStr testingStr filesep];
    if not(isfolder(warpBaseDir)); mkdir(warpBaseDir); end
    
    % TransformParameters of low-res registration to apply to high res data
    regiDir = [regiBaseDir ...
        'v6/non-rigid/testsGrad/Def_RegSeries_mask_LMall_gs16_bep1000_wLM0p1/'];
    lrRegistrationTransform = [regiDir 'TransformParameters.0.txt'];
    
    lr_warped_filename = [regiDir 'result.0.mha'];
end

% not needed, but can be useful to have
uct_regi_filename = [regiBaseDir 'volumes' filesep 'reco_b32.mha'];
atlas_regi_filename = [regiBaseDir 'volumes' filesep 'template_25.mha'];
atlas_seg_filename = [regiBaseDir 'volumes' filesep 'annotation_25.mha'];

% extended fixed mask used during registration
uct_mask_filename = [regiBaseDir 'volumes' filesep 'reco_b32_SE10.mha'];
atlas_mask_filename = [regiBaseDir 'volumes' filesep 'annotation_25_mask.mha'];

%% Defining ROIs to transform
% % - Open the fixed volume in ITK-SNAP
% % - Define ROI in terms of "Cursor positions (x,y,z)", which are pixels
% in image units [i.e. pixels]
% % - ROIs are given in form [x1,x2,y1,y2,z1,z2]

%% Transforming full image or series of ROIs 
% % Define a series of roi's to be transformed.
% % Each roi to be transformed needs the following in the form of cell arrays:
% % - a name (string in roiNameList)
% % - a roi vector ([x1,x2,y1,y2,z1,z2] in roiList)
% image units [i.e. pixels] in itk for fixed image

% fixed image is of size 1-456, 1-320, 1-528
if (1==0)
    % larger roi with full slice
    roiNameList = {'test_fullSliceroi'};
    roiList = {[1,456,1,320,65-5,65+5]};
    %roiNameList = {'test_fullSliceroi65','test_fullSliceroi75'};
    %roiList = {[1,456,1,320,60,70],[1,456,1,320,70,80]};
    %roiNameList = {'test_fullSliceroi5','test_fullSliceroi15'};
    %roiList = {[1,456,1,320,1,10],[1,456,1,320,10,20]};
    roiNameList = {'test_fullSliceroi15'};
    roiList = {[1,456,1,320,10,20]};
    %roiNameList = {'test_fullSliceroi5'};
    %roiList = {[1,456,1,320,1,10]};
    %roiNameList = {'test_fullSliceroi525'};
    %roiList = {[1,456,1,320,520,528]};
else
    % do full image, needs splitting
    roiNameList = {'full'};
    roiList = {[1,456,1,320,1,528]};
end

% name isotropic splitting different
if isotropicSplit
    isoStr = '_iso';
elseif optimizedSplit
    isoStr = '_opt';
else
    isoStr = '';
end
    
% output directories for tiffs
outBaseDir = [regiBaseDir 'tiffs' regStr isoStr testingStr filesep];
if ~exist(outBaseDir,'file')
    mkdir(outBaseDir)
end
outDir = [outBaseDir 'reco_b' num2str(res) '_warped' filesep];
if ~exist(outDir,'file')
    mkdir(outDir)
end

% output log file
dateStr = char(datetime);
dateStr(dateStr==':')='-';
dateStr(dateStr==' ')='-';

logFile = [warpBaseDir 'reco_b' num2str(res) '_' roiNameList{1}  isoStr '_log_' dateStr '.txt'];
diary(logFile)

%% output information what will be done
disp(['loadBaseDir    = ' loadBaseDir])
disp(['registration   = ' lrRegistrationTransform])
disp(['fixed image    = ' atlas_regi_filename])
disp(['moving image   = ' uct_regi_filename])
disp(['fixed segm.im. = ' atlas_regi_filename])
disp(' ')
disp(['fixed  resolution = ' num2str(ps_atlas_regi,'%.1f') ' um'])
disp(['moving resolution = ' num2str(ps_uct_regi,'%.1f') ' um ('  num2str(bf_regi) 'x)'])
disp(['out    resolution = ' num2str(ps_uct_1x*res,'%.1f') ' um ('  num2str(res) 'x)'])
disp(' ')
disp(['roiNameList{1} = ' roiNameList{1}])
disp(['roiList{1}     = ' num2str(roiList{1},'%d ')])
disp(['extendBy       = ' num2str(extendBy,'%d') ' pixels'])
disp(['MaxSizeGB      = ' num2str(MaxSizeGB,'%.2f') ' GB'])
disp(' ')
disp(['outDir warping = ' warpBaseDir])
disp(['outDir tiffs   = ' outDir])
disp(['logFile        = ' logFile])

%%
disp('press any key to continue')
disp(' ')
pause
