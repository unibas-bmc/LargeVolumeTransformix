%% Warping high-res with low-res transform
% % Pre-requisites:
% - uCT reconstructions were binned using ""
%   - bin 1 recos are tiff stacks in a directory "reco"
%   - other binnings are tiff stacks in directories "reco_bX"
% - registered uCT was produced as follows: 
%   - run  "BinningForRegistration.m" to go from b8 tiff to b32 tiff stack
%   - b32 tiff stack imported in FIJI/ImageJ, z-direction is flipped, and
%   exported as mha

% % Note about registration:
% The registraion currently used is a two-step registration: 1. manual
% pre-alignment, 2. automatic affine
% I give the TranformParameters files, but you will need to manually fix
% the paths for the affine TransformParameters:
% the field InitialTransformParametersFileName should give the full path to
% the manual prealignment parameters "ManualPrealignTParams.txt"
%
% If you would like to make your own manual pre-alignment, do the
% following:
% 1. open fixed and moving in ITK-SNAP
% 2. Tools -> Registration -> Manual
% 3. Manually align as best as possible using the interactive tool
% 4. Save transformation file (bottom, save icon)
% 5. Run a dummy elastix registration between moving and fixed (I provide
% dummy elastix parameters in
% ./example/registration/ManualPrealign/dummy_params.txt
% 6. Edit the resulting TransformParameters.txt, pasting the parameters
% from ITK-SNAP's saved file into the "TransformParameters" field and
% changing "CenterOfRotationPoint" to 0.0 0.0 0.0
% 7. You now have a prealignment TransformParameters file for running
% elastix with the -t0 flag

% % Naming conventions:
% - fixed image is "atlas"
% - moving image is "uct"

%% Toolboxes and functions
addpath('./utils/')

% % coordinate transforms
matlab2itk = @(coords) [coords(:,2)-1,coords(:,1)-1,coords(:,3)-1];
itk2matlab = @(coords) [coords(:,2)+1,coords(:,1)+1,coords(:,3)+1];

index_to_points = @(coords,ps) (coords-1)*ps;
points_to_index = @(coords,ps) round(coords/ps); % this is how transformix does it

lr_to_hr = @(coords,bf) (coords-1)*bf; % low resolution to high res coordinates

% % useful for outputting images
cropaboutcenter = @(im,roisize) im(floor(size(im,1)/2)-floor(roisize(2)/2)+1:floor(size(im,1)/2)+floor(roisize(2)/2),...
    floor(size(im,2)/2-roisize(1)/2)+1:floor(size(im,2)/2+roisize(1)/2));
cropto = @(im,roi) im(roi(3):roi(4),roi(1):roi(2));
cropto3d = @(im,roi) im(roi(3):roi(4),roi(1):roi(2),roi(5):roi(6));
croptostack = @(im,roi) im(roi(3):roi(4),roi(1):roi(2),:);
prepimage = @(im,vr) uint8(255*(double(im)-vr(1))/(vr(2)-vr(1)));
prepimagei16 = @(im,vr) uint16((2^16-1)*(double(im)-vr(1))/(vr(2)-vr(1)));

%% Directory structure
% base directory where high resolution tiff stacks are stored
loadBaseDir = '/media/griffin/anatomix_mousebrain/mouse4_perf_eth/';
% for bin 1, tiffs should be here: [loadBaseDir 'reco' filesep];
% for any other binning, tiffs should be here: [loadBaseDir 'reco_bX' filesep];

% base directory where warped images will be stored
warpBaseDir = './example/registration/warping/';
if not(isfolder(warpBaseDir)); mkdir(warpBaseDir); end

% location of registration files (used only in this section)
regiBaseDir = './example/registration/';

% TransformParameters of low-res registration to apply to high res data
lrRegistrationTransform = [regiBaseDir ...
    'affine0/TransformParameters.0.txt'];

% not needed, but can be useful to have
uct_regi_filename = [regiBaseDir 'volumes/reco_b32.mha'];
atlas_regi_filename = [regiBaseDir 'volumes/template_25.mha'];

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

%% Defining ROIs to transform
% % - Open the fixed volume in ITK-SNAP
% % - Define ROI in terms of "Cursor positions (x,y,z)", which are pixels
% in image units [i.e. pixels]
% % - ROIs are given in form [x1,x2,y1,y2,z1,z2]

%% Transforming small ROI (i.e. no stitching necessary)
% % Define a series of roi's to be transformed. 
% % Each roi to be transformed needs the following in the form of cell arrays:
% % - a name (string in roiNameList)
% % - a roi vector ([x1,x2,y1,y2,z1,z2] in roiList)
% % - a list of resolutions to be transformed ([r1,r2,...] in resLists)
roiNameList = {'test_smallroi'};
% image units [i.e. pixels] in itk for fixed image
roiList = {[100,180,60,140,65-5,65+5]}; 
% which resolutions you want (i.e. binning factor) -- currently works for 1, 2, 4, and 8
resLists = {[8]}; 

% Do you want to remove the unwarped high res volumes used as an
% intermediate to get the warped result? (useful to have for debugging)
removeUnwarped = 0;

% % Loop over ROIs to tranform
for r = 1:length(roiNameList)
    roi = roiList{r};
    roiName = roiNameList{r};
    resList = resLists{r};
    
    %% Transform points from registered to inital frames
    % directory to put all point transformations
    roiBaseDir = [warpBaseDir roiName filesep];
    if not(isfolder(roiBaseDir)); mkdir(roiBaseDir); end

    % transform coordinates to bounding box then from matlab to itk
    iniroi_itk = roi; % initial region of interest as defined in itk (image units)
    iniroi_bb_itk = RoiToBoundaryCoordinates(iniroi_itk); % get bounding box
    iniroi_bb_itk_points = index_to_points(iniroi_bb_itk,ps_atlas_regi/ps_fac);
    
    % write these to a file which transformix can use
    fname = [roiBaseDir 'iniroi_bbpoints.txt'];
    writePointsFile(fname,iniroi_bb_itk_points,'point');
    
    % run transformix on this with the forward transform
    calltransformix(fname,lrRegistrationTransform,roiBaseDir);
    
    % read the transformed points
    readName = [roiBaseDir 'outputpoints.txt'];
    outPoints = readTransformedPointsFile(readName);
    traroi_bb1_itk_points = outPoints.OutputPoint;
    traroi_bb1_itk = outPoints.OutputIndexFixed;
    
    %% Loop over desired resolutions, load volumes, and transform
    for i = 1:length(resList)
        % which resolution?
        res = resList(i);
        
        % directory for this resolution
        thisDir = [roiBaseDir num2str(res) 'x' filesep];
        if not(isfolder(thisDir)); mkdir(thisDir); end
        
        % directory where reco tiff stacks live
        if res == 1
            thisLoadDir = [loadBaseDir 'reco' filesep];
        else
            thisLoadDir = [loadBaseDir 'reco_b' num2str(res) filesep];
        end
        
        % size of reco at this resolution
        [sx,sy,sz] = GetStackSize(thisLoadDir);
        
        % ROI transformed into points in moving image
        roi_uct_regires = points_to_index(traroi_bb1_itk_points,ps_uct_regi/ps_fac);
        roi_uct_thisres = points_to_index(traroi_bb1_itk_points,(ps_uct_1x*res)/ps_fac);
        
        % get bounding box for this skewed roi
        loadroi_uct_thisres = CoordsToBoundingboxRoi(roi_uct_thisres);
        % flip z because tiff stack is flipped z w.r.t. mha
        loadroi_uct_thisres(5:6) = sz-[loadroi_uct_thisres(6),loadroi_uct_thisres(5)];
        
        % size of roi
        nvox = abs(diff(loadroi_uct_thisres));
        nvox = prod(nvox([1,3,5])); % number of voxels in unwarped roi
        GB = nvox*2e-9; % size in GB of unwarped roi
        fprintf('Loading unwarped roi for ROI:%s, res:%d\nUnwarped volume will be %f GB in size\n',...
            roiName,res,GB)
        
        vol_uw = flipud(rot90(flip(stackreader(thisLoadDir,loadroi_uct_thisres),3)));
        
        % writing unwarped volume to mha
        fname = sprintf('%s%sb%d_unwarped.mha',thisDir,roiName,res);
        resolution = (ps_uct_1x*res)*ones(1,3)/ps_fac;
        offset = min(traroi_bb1_itk_points,[],1);
        writeMetaImageFile(fname,vol_uw,resolution,offset);
        
        % get target coordinates in hr registered frame
        tmp1 = min(iniroi_bb_itk_points,[],1);
        tmp2 = max(iniroi_bb_itk_points,[],1);
        xax = tmp1(1):(ps_uct_1x*res/ps_fac):tmp2(1);
        yax = tmp1(2):(ps_uct_1x*res/ps_fac):tmp2(2);
        zax = tmp1(3):(ps_uct_1x*res/ps_fac):tmp2(3);
        outSize = [length(xax),length(yax),length(zax)];
        outRes = (ps_uct_1x*res/ps_fac)*[1,1,1];
        outOrig = min(iniroi_bb_itk_points,[],1);
        
        % generate transform file for transformix
        thisTransform = sprintf('%s%s_b%d_TransformParameters.txt',...
            thisDir,roiName,res);
        
        lineno1 = 12; % image size line
        strnew1 = ['(Size ' num2str(outSize(1)) ' ' num2str(outSize(2)) ' ' num2str(outSize(3)) ')'];
        lineno2 = 14; % image resolution line
        strnew2 = ['(Spacing ' num2str(outRes(1)) ' ' num2str(outRes(2)) ' ' num2str(outRes(3)) ')'];
        lineno3 = 15; % origin line
        strnew3 = ['(Origin ' num2str(outOrig(1)) ' ' num2str(outOrig(2)) ' ' num2str(outOrig(3)) ')'];
        tparamreplaceline(lrRegistrationTransform,thisTransform,...
            lineno1,strnew1,lineno2,strnew2,lineno3,strnew3);
        
        % call transformix to warp volume
        calltransformix(fname,thisTransform,thisDir);
        
        % (optionally) remove unwarped volume
        if removeUnwarped == 1
            delete(fname);
        end
    end
end

%% Transforming large ROI (i.e. stitching necessary)
% % Bigger to-do:
% - run transformix on gpu
% - can we load faster if we tiling the tiffs? see tileData = readEncodedTile(t,tileNumber)
% % Smaller to-do:
% - padded stack reader
% - padding sub-roi for non-rigid warping --- pad either the initial roi 
% (only for loading, not for target volume) or pad the loadroi
%    -- be sure to account for different origin and sizes of padded rois
%       when saving unwarped roi

MaxSizeGB = 0.02; % [GB] break into rois no larger than MaxSizeGB 
% be conservative, this won't be exact and is calculated in the registered
% frame. That means in the case a skewed roi in the floating frame will
% lead to bigger sizes. Also, some overhead should be given for 

roiNameList = {'test1'};
%roiList = {[1,456,1,320,266-20,266+20]}; % image units [i.e. pixels] in itk for fixed image
roiList = {[150,170,80,100,220,240]}; % image units [i.e. pixels] in itk for fixed image
resLists = {[2]}; % which resolutions you want (i.e. binning factor) -- currently works for 1, 2, 4, and 8

removeUnwarped = 1; % 1: deletes unwarped after warping, 0: leaves unwarped on disk

for r = 1:length(roiNameList)
    % which roi?
    roi = roiList{r};
    roiName = roiNameList{r};
    resList = resLists{r};
    
    % size of roi
    tmp = diff(roi);
    roiSizeAtlas = tmp([1,3,5]);
        
    % base directory: where to put this transformed roi
    roiBaseDir = [warpBaseDir roiName filesep];
    if not(isfolder(roiBaseDir)); mkdir(roiBaseDir); end

    %% Loop over desired resolutions
    for i = 1:length(resList)
        % which resolution?
        res = resList(i);
        
        % directory for this resolution
        thisDir = [roiBaseDir num2str(res) 'x' filesep];
        if not(isfolder(thisDir)); mkdir(thisDir); end
        
        % directory where reco tiff stacks live
        if res == 1
            thisLoadDir = [loadBaseDir 'reco' filesep];
        else
            thisLoadDir = [loadBaseDir 'reco_b' num2str(res) filesep];
        end
        
        % size of reco at this resolution
        [sx,sy,sz] = GetStackSize(thisLoadDir);
        
        % size of roi at this resolution
        roiSizeThisRes = ceil(roiSizeAtlas*(ps_atlas_25/res));
        roiSizeGB = 2*prod(roiSizeAtlas*(ps_atlas_25/res))/1e9;
        
        % will break into grid of subvolumes: gridFac x gridFac x gridFac
        gridFac = ceil((roiSizeGB/MaxSizeGB)^(1/3));
        %gridFac = nextpow2((roiSizeGB/MaxSizeGB)^(1/3));
        [newRois,roiIndices] = SubRoisFromRoi(roi,gridFac);
        % % Note: may want to come up with scheme for case where non-cubic roi
        
        %% Loop over sub-rois
        for r2 = 1:length(newRois)
            thisroi = newRois{r2};
            
            % sub-roi directories
            tmp = roiIndices{r2};
            thisSubROIDir = [sprintf('%s%d_%d_%d',thisDir,tmp(1),tmp(2),tmp(3)) filesep];
            if not(isfolder(thisSubROIDir)); mkdir(thisSubROIDir); end
            
            %% Transform points from registered to inital frames
            % directory to put all point transformations
            pointsfileLocation = thisSubROIDir;
            
            % transform coordinates to bounding box then from matlab to itk
            iniroi_itk = thisroi; % initial region of interest as defined in itk (image units)
            iniroi_bb_itk = RoiToBoundaryCoordinates(iniroi_itk); % get bounding box
            iniroi_bb_itk_points = index_to_points(iniroi_bb_itk,ps_atlas_regi/ps_fac);
            
            % write these to a file which transformix can use
            fname = [thisSubROIDir 'iniroi_bbpoints.txt'];
            writePointsFile(fname,iniroi_bb_itk_points,'point');
            
            % run transformix on this with the forward transform
            calltransformix(fname,lrRegistrationTransform,pointsfileLocation);
            
            % read the transformed points
            readName = [pointsfileLocation 'outputpoints.txt'];
            outPoints = readTransformedPointsFile(readName);
            traroi_bb1_itk_points = outPoints.OutputPoint;
            traroi_bb1_itk = outPoints.OutputIndexFixed;
            
            
            roi_uct_regires = points_to_index(traroi_bb1_itk_points,ps_uct_regi/ps_fac);
            roi_uct_thisres = points_to_index(traroi_bb1_itk_points,(ps_uct_1x*res/ps_fac));
            
            loadroi_uct_thisres = CoordsToBoundingboxRoi(roi_uct_thisres);
            loadroi_uct_thisres(5:6) = sz-[loadroi_uct_thisres(6),loadroi_uct_thisres(5)];
            
            nvox = abs(diff(loadroi_uct_thisres));
            nvox = prod(nvox([1,3,5])); % number of voxels in unwarped roi
            GB = nvox*2e-9; % size in GB of unwarped roi
            
            vol_uw = flipud(rot90(flip(stackreader(thisLoadDir,loadroi_uct_thisres),3)));
            
            % writing unwarped volume to mha
            fname = [sprintf('%s%s_%dx_%d_%d_%d',thisSubROIDir,roiName,res,tmp(1),tmp(2),tmp(3)) '_unwarped.mha'];
            resolution = (ps_uct_1x*res)*ones(1,3)/ps_fac;
            offset = min(traroi_bb1_itk_points,[],1);
            writeMetaImageFile(fname,vol_uw,resolution,offset);
            
            % get target coordinates in hr registered frame
            tmp1 = min(iniroi_bb_itk_points,[],1);
            tmp2 = max(iniroi_bb_itk_points,[],1);
            xax = tmp1(1):(ps_uct_1x*res/ps_fac):tmp2(1);
            yax = tmp1(2):(ps_uct_1x*res/ps_fac):tmp2(2);
            zax = tmp1(3):(ps_uct_1x*res/ps_fac):tmp2(3);
            outSize = [length(xax),length(yax),length(zax)];
            outRes = (ps_uct_1x*res)*[1,1,1]/ps_fac;
            outOrig = min(iniroi_bb_itk_points,[],1);
            
            % generate transform file for transformix
            thisTranform = sprintf('%s%s_%dx_%d_%d_%d_TransformParameters.txt',thisSubROIDir,roiName,res,tmp(1),tmp(2),tmp(3));
            
            lineno1 = 12; % image size line
            strnew1 = ['(Size ' num2str(outSize(1)) ' ' num2str(outSize(2)) ' ' num2str(outSize(3)) ')'];
            lineno2 = 14; % image resolution line
            strnew2 = ['(Spacing ' num2str(outRes(1)) ' ' num2str(outRes(2)) ' ' num2str(outRes(3)) ')'];
            lineno3 = 15; % origin line
            strnew3 = ['(Origin ' num2str(outOrig(1)) ' ' num2str(outOrig(2)) ' ' num2str(outOrig(3)) ')'];
            tparamreplaceline(lrRegistrationTransform,thisTranform,...
                lineno1,strnew1,lineno2,strnew2,lineno3,strnew3);
            
            % call transformix to warp volume
            calltransformix(fname,thisTranform,thisSubROIDir);
            
            % (optionally) remove unwarped volume
            if removeUnwarped == 1
                delete(fname);
            end
        end
        
        %% Load warped roi mha files and write out full tiffs
        for bz = 1:gridFac
            
            
            
            
            
            
            
        end
        
        %% (if desired) delete the warped volume files to free up disk space
                
    end
end

