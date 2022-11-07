%% Warping high-res with low-res transform
% % Pre-requisites:
% - uCT reconstructions were binned using ""
%   - bin 1 recos are tiff stacks in a directory "reco"
%   - other binnings are tiff stacks in directories "reco_bX"
% - registered uCT was produced as follows:
%   - run  "BinningForRegistration.m" to go from b8 tiff to b32 tiff stack
%   - b32 tiff stack imported in FIJI/ImageJ, z-direction is flipped
%     (image->transform->flip Z), voxel size set (image->properties) and
%     saved as mha

% % Note about registration:
% The registration currently used is a two-step registration: 1. manual
% pre-alignment, 2. automatic affine
% The TranformParameters files are provided, with the affine TransformParameters
% file including the relative path to the manual transformation for the
% field InitialTransformParametersFileName, i.e.
% "./example/registration/ManualPrealign/ManualPrealignTParams.txt"
% This needs adjusting if directory structure changes
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

% to-do:
% - run transformix on gpu
% - can we load faster if we tiling the tiffs? see tileData = readEncodedTile(t,tileNumber)
% - logic of setting origin for all resolutions such that data aligns. This depends on binning and
%   interpretation of next processing step (e.g. itksnap vs. elastix vs. neuroglancer)
% - write-out data in precomputed neuroglancer format

% load data path and parameter definitions
DataParameterDefinition

%% Toolboxes and functions
addpath('./utils/')

% record start time
startTime = tic;

% % coordinate transforms
matlab2itk = @(coords) [coords(:,2)-1,coords(:,1)-1,coords(:,3)-1];
itk2matlab = @(coords) [coords(:,2)+1,coords(:,1)+1,coords(:,3)+1];

if (1==1)
    % how transformix does it
    index_to_points = @(i,ps) i*ps;         % index i to world x coords 
    points_to_index = @(x,ps) round(x/ps);  % world x to index i coords
else
    % this is how ITK does it, index i=[1 1 1] is x=[0 0 0] - more differences
    index_to_points = @(i,ps) (i-1)*ps;        % index i to world x coords
    points_to_index = @(x,ps) round((x/ps)+1); % world x to index i coords
end

lr_to_hr = @(coords,bf) (coords-1)*bf; % low resolution to high res coordinates

% % useful for outputting images
cropaboutcenter = @(im,roisize) im(floor(size(im,1)/2)-floor(roisize(2)/2)+1:floor(size(im,1)/2)+floor(roisize(2)/2),...
    floor(size(im,2)/2-roisize(1)/2)+1:floor(size(im,2)/2+roisize(1)/2));
cropto = @(im,roi) im(roi(3):roi(4),roi(1):roi(2));
cropto3d = @(im,roi) im(roi(3):roi(4),roi(1):roi(2),roi(5):roi(6));
croptostack = @(im,roi) im(roi(3):roi(4),roi(1):roi(2),:);
prepimage = @(im,vr) uint8(255*(double(im)-vr(1))/(vr(2)-vr(1)));
prepimagei16 = @(im,vr) uint16((2^16-1)*(double(im)-vr(1))/(vr(2)-vr(1)));

if doTesting
    % for TESTING - warp low-res (b32) moving image
    lr_warped_filename = [regiBaseDir 'volumes/result.mha'];
    if ~exist(lr_warped_filename,'file')
        calltransformix(uct_regi_filename,lrRegistrationTransform,[regiBaseDir 'volumes/']);
    end
    warped_full_lr = mha_read_volume(lr_warped_filename);

    % Do you want to remove the unwarped high res volumes used as an
    % intermediate to get the warped result? (useful to have for debugging)
    removeUnwarped = 0;
else
    removeUnwarped = 1;
end

% % Loop over ROIs to tranform
for r = 1:length(roiNameList)
    roi = roiList{r};
    roiName = roiNameList{r};
    
    %% Transform points from registered to inital frames
    % directory to put all point transformations
    roiBaseDir = [warpBaseDir roiName filesep];
    if not(isfolder(roiBaseDir)); mkdir(roiBaseDir); end
    
    % transform coordinates to bounding box then from matlab to itk
    iniroi_itk = roi; % initial region of interest as defined in itk (image units)
    iniroi_bb_itk = RoiToBoundaryCoordinates(iniroi_itk); % get bounding box
    
    %% load volumes, and transform
    
    % express bounding box etc. w.r.t. HR pixel locations
    % world coordinates of LR fixed image
    iniroi_bb_itk_LR_points = index_to_points(iniroi_bb_itk,ps_atlas_regi/ps_fac);
    % convert these to indices of HR image
    % true inverse and rounding
    iniroi_bb_itk_HR = points_to_index(iniroi_bb_itk_LR_points,(ps_uct_1x*res)/ps_fac);
    % convert back to points in LR image
    iniroi_bb_itk_points = index_to_points(iniroi_bb_itk_HR,(ps_uct_1x*res)/ps_fac);
    
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
    
    % directory for this resolution
    thisDir = [roiBaseDir num2str(res) 'x' filesep];
    if not(isfolder(thisDir)); mkdir(thisDir); end
    
    % directory where reco tiff stacks live
    if res == 1
        thisLoadDir = [loadBaseDir 'reco' filesep];
    else
        thisLoadDir = [loadBaseDir 'reco_b' num2str(res) filesep];
    end
    % extract only once list of files as this costs time
    fDir=dir([thisLoadDir '*.tif']);
    
    % size of reco at this resolution
    [sx,sy,sz] = GetStackSize(fDir);
    
    % ROI transformed into points in moving image
    roi_uct_regires = points_to_index(traroi_bb1_itk_points,ps_uct_regi/ps_fac);
    roi_uct_thisres = points_to_index(traroi_bb1_itk_points,(ps_uct_1x*res)/ps_fac);
    
    % get bounding box for this skewed roi
    loadroi_uct_thisres = CoordsToBoundingboxRoi(roi_uct_thisres);
    % flip z because tiff stack is flipped z w.r.t. mha
    loadroi_uct_thisres(5:6) = sz-[loadroi_uct_thisres(6),loadroi_uct_thisres(5)];
   
    % increase roi of unwarped image to avoid interpolation
    % artifacts
    loadroi_uct_thisres = loadroi_uct_thisres + extendBy*[-1 1 -1 1 -1 1];
    extOffset = extendBy*[-1 -1 -1];
    
    % size of roi
    nvox = abs(diff(loadroi_uct_thisres));
    nvox = prod(nvox([1,3,5])+[1 1 1]); % number of voxels in unwarped roi
    GB = nvox*2e-9; % size in GB of unwarped roi type int16, factor 2 checked
    fprintf('Unwarped volume for ROI:%s, res:%d will be %.2f GB in size\n',roiName,res,GB)
   
    if GB<MaxSizeGB || doTesting
        % can do full region, do not need to split it
        % only produce this for testing
        
        % int16 for b8, while uint16 for b32!
        vol_uw = flipud(rot90(flip(stackreader(fDir,loadroi_uct_thisres),3)));
        if isa(vol_uw,'uint16')
            % make it int16
            newMin=double(intmin('int16'));
            % write out as int16
            vol_uw=int16(double(vol_uw)+newMin);
        end
        
        if showROIBox && doTesting
            % visualize region to extract
            nvoxV = abs(diff(loadroi_uct_thisres));
            BB = [loadroi_uct_thisres([1 3 5]) nvoxV([1,3,5])];  % bounding box
            centerBB = BB(1:3)+0.5*BB(4:6);
            figure(77)
            fh=drawCuboid(BB(4:6)',centerBB',[0;0;0],'g');
            hold on
            plot3(roi_uct_thisres(:,1),roi_uct_thisres(:,2),sz-roi_uct_thisres(:,3),'ro')
            %plot3(roi_uct_thisres(1,1),roi_uct_thisres(1,2),sz-roi_uct_thisres(1,3),'r*')
            %plot3(centerBB(1,1),centerBB(1,2),centerBB(1,3),'ko')
            hold off
            view(3)
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end
        
        % writing unwarped volume to mha
        fnameUW = sprintf('%s%sb%d_unwarped.mha',thisDir,roiName,res);
        outRes = (ps_uct_1x*res)*ones(1,3)/ps_fac;
        % this offset makes the unwarped volumes aligned in itksnap for b32
        offsetForITK = [-outRes(1:2) 0];
        % derive offset directly from index from roi
        offset = index_to_points(min(roi_uct_thisres)+extOffset,(ps_uct_1x*res)/ps_fac) + offsetForITK;
        writeMetaImageFile(fnameUW,vol_uw,outRes,offset);
        
        % get target coordinates in hr registered frame
        tmp1 = min(iniroi_bb_itk_points,[],1);
        tmp2 = max(iniroi_bb_itk_points,[],1);
        xax = tmp1(1):outRes(1):tmp2(1)-outRes(1);  % pixel coordinates
        yax = tmp1(2):outRes(2):tmp2(2)-outRes(2);
        zax = tmp1(3):outRes(3):tmp2(3)-outRes(3);
        outSize = [length(xax),length(yax),length(zax)];
        outOrig = tmp1-outRes; % itk origin, needs subtracting of outRes, ok b32
        
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
        calltransformix(fnameUW,thisTransform,thisDir);
        
        % (optionally) remove unwarped volume
        if removeUnwarped == 1
            delete(fnameUW);
        else
            % check that unwarped images align
            disp(['itksnap -o ' regiBaseDir 'volumes/reco_b32.mha -g ' fnameUW ' &'])
            
            disp(['roiUW loadroi R0 ' num2str(loadroi_uct_thisres,'%d ')])
            disp(['roiUW minroi  R0 ' num2str(min(roi_uct_thisres)+extOffset,'%d ')])
            disp(['roiUW offset  R0 ' num2str(offset,'%.3f ')])
            
            minIdx = points_to_index( tmp1,(ps_uct_1x*res/ps_fac));
            maxIdx = points_to_index( tmp2,(ps_uct_1x*res/ps_fac));
            disp(['minCoo R0 ' num2str(tmp1,'%.2f ')])
            disp(['maxCoo R0 ' num2str(tmp2,'%.2f ')])
            disp(['minIdx R0 ' num2str(minIdx,'%d ')])
            disp(['maxIdx R0 ' num2str(maxIdx,'%d ')])
            %keyboard
        end
        
        if doTesting && r==length(roiNameList)
            %% visually compare results
            disp(['x 0 ' num2str(roi(1:2),'%.1f ') 'px, ' num2str([xax(1) xax(end)],'%.3f ') 'um HR, res '  num2str(outRes(1)) 'um, size ' num2str(outSize(1)) 'px'])
            
            zTest=1;   % show selected slice
            
            % show warped HR roi
            fname = [thisDir 'result.mha'];
            warpSmallRoiMHA = mha_read_volume(fname);
            figure(10)
            imagesc(warpSmallRoiMHA(:,:,zTest))
            title(['warped HR, res. ' num2str(res*ps_uct_1x) 'um'])
            axis equal tight
            colormap('gray')
            
            % downsample to size of fixed image resolution
            % NOT taking care of different voxel positions!
            figure(12)
            tmpI=imresize3(warpSmallRoiMHA,res*ps_uct_1x/ps_atlas_regi);
            imagesc(tmpI(:,:,zTest))
            title(['downsampled warped HR, res. ' num2str(ps_atlas_regi) 'um'])
            axis equal tight
            colorbar
            colormap('gray')
            
            % same roi of warped LR
            figure(11)
            tmpLR=warped_full_lr(roi(1):roi(2)-1,roi(3):roi(4)-1,roi(5):roi(6)-1);
            imagesc(tmpLR(:,:,zTest))
            title(['warped LR, res. ' num2str(ps_atlas_regi) 'um'])
            axis equal tight
            colorbar
            colormap('gray')
            
            binV=linspace(double(intmin('int16')),double(intmax('int16')),100);
            %binV=linspace(double(min([tmpI(:); tmpLR(:)])),double(max([tmpI(:); tmpLR(:)])),100);
            figure(18)
            subplot(4,1,1)
            histogram(double(vol_uw(:)),binV);
            title(['unwarped HR, res. ' num2str(res*ps_uct_1x) 'um'])
            subplot(4,1,2)
            histogram(double(warpSmallRoiMHA(:)),binV);
            title(['warped HR, res. ' num2str(res*ps_uct_1x) 'um'])
            subplot(4,1,3)
            histogram(double(tmpI(:)),binV);
            title(['downsampled warped HR, res. ' num2str(ps_atlas_regi) 'um'])
            subplot(4,1,4)
            histogram(double(tmpLR(:)),binV);
            title(['warped LR, res. ' num2str(ps_atlas_regi) 'um'])
            
        end
        disp(['itksnap -o ' regiBaseDir 'volumes/result.mha -g ' regiBaseDir 'warping/' roiName '/' num2str(res) 'x/result.mha &'])
        % small enough to do it in one go
    end
    % for all regions
end

%% Transforming large ROI (i.e. stitching necessary)

if doTesting
    % use same last roi as before to compare
    removeUnwarped = 0; % 1: deletes unwarped after warping, 0: leaves unwarped on disk
else
    removeUnwarped = 1; % 1: deletes unwarped after warping, 0: leaves unwarped on disk
end
for r=1:length(roiNameList)
    roiNameList{r} = [roiNameList{r} '_splitted' isoStr];
end
% be conservative, this won't be exact and is calculated in the registered
% frame. That means in the case a skewed roi in the floating frame will
% lead to bigger sizes. Also, some overhead should be given for

for r = 1:length(roiNameList)
    % which roi?
    roi = roiList{r};
    roiName = roiNameList{r};
    
    % size of roi
    tmp = diff(roi);
    roiSizeAtlas = tmp([1,3,5]);
    
    % base directory: where to put this transformed roi
    roiBaseDir = [warpBaseDir roiName filesep];
    if not(isfolder(roiBaseDir)); mkdir(roiBaseDir); end
    
    % directory for this resolution
    thisDir = [roiBaseDir num2str(res) 'x' filesep];
    if not(isfolder(thisDir)); mkdir(thisDir); end
    
    % directory where reco tiff stacks live
    if res == 1
        thisLoadDir = [loadBaseDir 'reco' filesep];
    else
        thisLoadDir = [loadBaseDir 'reco_b' num2str(res) filesep];
    end
    % extract only once list of files as this costs time
    fDir=dir([thisLoadDir '*.tif']);
    
    % size of reco at this resolution
    [sx,sy,sz] = GetStackSize(fDir);
    
    % size of whole roi at HR
    roiSizeThisRes = ceil(roiSizeAtlas*(ps_atlas_25/res));
    % estimate size of warped HR image
    roiSizeGB = prod(roiSizeThisRes)*8e-9;
    fprintf('Warped volume for whole ROI:%s, res:%d will be %.2f GB in size\n',roiName,res,roiSizeGB)
   
    % will break into grid of subvolumes: gridFacV(1) x gridFacV(2) x gridFacV(3)
    if isotropicSplit
        % isotropic splitting, i.e. same number of grids
        % taking longer due to reading many mha's for each output tiff slice
        gridFac = ceil((roiSizeGB/MaxSizeGB)^(1/3));
        gridFacV = [gridFac gridFac gridFac];
    else
        % ideally keep full in-plane slices
        gridFacXY=1;
        gridFacZ=ceil(roiSizeGB/(gridFacXY*gridFacXY*MaxSizeGB)); 
        gridFacV = [gridFacXY gridFacXY gridFacZ];
    end
    disp(['splitting into ' num2str(gridFacV,'%d ') ' subimages'])
    
    % keep exact pixel locations as for roi
    % rather split up pixel coordinate vector
    % determine coordinates for roi, 8 points of bounding box
    roi_bb_itk = RoiToBoundaryCoordinates(roi); % get bounding box
    % world coordinates of LR fixed image
    roi_bb_itk_LR_points = index_to_points(roi_bb_itk,ps_atlas_regi/ps_fac);
    % convert these to indices of HR image
    roi_bb_itk_HR = points_to_index(roi_bb_itk_LR_points,(ps_uct_1x*res)/ps_fac);
    roiHR([1 3 5])=min(roi_bb_itk_HR);
    roiHR([2 4 6])=max(roi_bb_itk_HR);
    % split high-resolution coordinates, all integer
    [newRois,roiIndices] = SubRoisFromRoiHR(roiHR,gridFacV);
    
    % % Note: may want to come up with scheme for case where non-cubic roi
    noRois = length(newRois);
    
    % save coordinates
    minCoord = zeros(3,noRois);
    maxCoord = zeros(3,noRois);
    
    %% Loop over sub-rois
    for r2 = 1:noRois
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
        iniroi_bb_itk_points = index_to_points(iniroi_bb_itk,(ps_uct_1x*res)/ps_fac);
        
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
        
        % increase roi of unwarped image to avoid interpolation artifacts
        loadroi_uct_thisres = loadroi_uct_thisres + extendBy*[-1 1 -1 1 -1 1];
        extOffset = extendBy*[-1 -1 -1];
        
        if r2==1
            % inform about file sizes
            nvox = abs(diff(loadroi_uct_thisres));
            nvox = prod(nvox([1,3,5])+[1 1 1]); % number of voxels in unwarped roi
            GB = nvox*2e-9; % size in GB of unwarped roi type int16, factor 2 checked
            fprintf('Unwarped volume for ROI:%s %d_%d_%d, res:%d will be %.2f GB in size\n',roiName,tmp,res,GB)
        
            % estimate size of warped HR roi image
            nvox = abs(diff(thisroi));
            nvox = prod(nvox([1,3,5])); % number of voxels in warped HR roi
            splitRoiSizeGB = nvox*2e-9; % size in GB of warped roi type int16, factor 2 checked
            fprintf('Warped volume for ROI:%s %d_%d_%d, res:%d will be %.2f GB in size\n',roiName,tmp,res,splitRoiSizeGB)
        end
   
        % extract unwarped moving image region
        vol_uw = flipud(rot90(flip(stackreader(fDir,loadroi_uct_thisres),3)));
        if isa(vol_uw,'uint16')
            % make it int16
            newMin=double(intmin('int16'));
            vol_uw=int16(double(vol_uw)+newMin);
        end
        
        if doTesting && r2==1 && showROIBox
            % visualize unwarped region to extract
            nvoxV = abs(diff(loadroi_uct_thisres));
            BB = [loadroi_uct_thisres([1 3 5]) nvoxV([1,3,5])];  % bounding box
            centerBB = BB(1:3)+0.5*BB(4:6);
            figure(77)
            hold on
            drawCuboid(BB(4:6)',centerBB',[0;0;0],'b');
            hold on
            %plot3(roi_uct_thisres(1,1),roi_uct_thisres(1,2),sz-roi_uct_thisres(1,3),'mo')
            plot3(roi_uct_thisres(:,1),roi_uct_thisres(:,2),sz-roi_uct_thisres(:,3),'k*')
            view(3)
            xlabel('x')
            ylabel('y')
            zlabel('z')
            axis equal
            legend('fullROI','fullPoints','1 1 1','Points 1 1 1')
            title('ROIs in unwarped image')
        end
        
        % writing unwarped volume to mha
        fnameUW = [sprintf('%s%s_%dx_%d_%d_%d',thisSubROIDir,roiName,res,tmp(1),tmp(2),tmp(3)) '_unwarped.mha'];
        outRes = (ps_uct_1x*res)*ones(1,3)/ps_fac;
        % this offset makes the unwarped volumes aligned for b32
        offsetForITK = (ps_uct_1x*res)/ps_fac*[-1 -1 0];
        % derive offset directly from index from roi
        offset = index_to_points(min(roi_uct_thisres)+extOffset,(ps_uct_1x*res)/ps_fac) + offsetForITK;
        writeMetaImageFile(fnameUW,vol_uw,outRes,offset);
        
        % get target coordinates in hr registered frame
        tmp1 = min(iniroi_bb_itk_points,[],1);
        tmp2 = max(iniroi_bb_itk_points,[],1);
        xax = tmp1(1):outRes(1):tmp2(1)-outRes(1);
        yax = tmp1(2):outRes(2):tmp2(2)-outRes(2);
        zax = tmp1(3):outRes(3):tmp2(3)-outRes(3);
        outSize = [length(xax),length(yax),length(zax)];
        outOrig = tmp1-outRes;   % itk origin, needs subtracting of outRes, ok for b32
        
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
        calltransformix(fnameUW,thisTranform,thisSubROIDir);
        
        % save target coordinates in hr registered frame
        minCoord(:,r2)=tmp1;
        maxCoord(:,r2)=tmp2-outRes;
        
        % (optionally) remove unwarped volume
        if removeUnwarped == 1
            delete(fnameUW);
        else
            if r2==1 %round(gridFac*gridFac/2)
                disp(['itksnap -o ' regiBaseDir 'volumes/reco_b32.mha -g ' fnameUW ' &'])
                
                disp(['roiUW loadroi R' num2str(r2) ' ' num2str(loadroi_uct_thisres,'%d ')])
                disp(['roiUW minroi  R' num2str(r2) ' ' num2str(min(roi_uct_thisres)+extOffset,'%d ')])
                disp(['roiUW offset  R' num2str(r2) ' ' num2str(offset,'%.3f ')])
                
                minIdx = points_to_index( tmp1,(ps_uct_1x*res/ps_fac));
                maxIdx = points_to_index( tmp2,(ps_uct_1x*res/ps_fac));
                disp(['minCoo R' num2str(r2) ' ' num2str(tmp1,'%.2f ')])
                disp(['maxCoo R' num2str(r2) ' ' num2str(tmp2,'%.2f ')])
                disp(['minIdx R' num2str(r2) ' ' num2str(minIdx,'%d ')])
                disp(['maxIdx R' num2str(r2) ' ' num2str(maxIdx,'%d ')])
                disp(strnew3)
                %keyboard
            end
        end
        
    end
    
    %% Load warped roi mha files and write out full tiff slices
    % ideally want rois to be large in-plane and small out of plane
    % to minimize reading overhead
    %
    % roi: fixed image ROI BB [in vxs] to be filled, but at
    % higher resolution
    %
    % newRois{}: cells with fixed image sub-roi BB [in vxs]
    %
    % TESTING with full 3D patches
    part1Time = toc(startTime);
    startTiffTime=tic;
     
    % loop over grid divisions in z
    for gridZ = 1:gridFacV(3)
        gridZoffset = (gridZ-1)*gridFacV(1)*gridFacV(2);
        
        % determine HR pixel coordinates
        minIdx = points_to_index( minCoord(:,1+gridZoffset),(ps_uct_1x*res/ps_fac));
        maxIdx = points_to_index( maxCoord(:,1+gridZoffset),(ps_uct_1x*res/ps_fac));
        zV = minIdx(3):maxIdx(3);
        
        disp(['gridZ ' num2str(gridZ) ', creating ' num2str(length(zV)) ' tiffs'])
   
        % loop over all z slices per gridZ
        for zIdx=1:length(zV)
            % loop over in-plane grid division
            for r1 = 1:min(gridFacV(1)*gridFacV(2),noRois)
                %
                r2=r1+gridZoffset;
                % voxel coordinates of ROI
                thisroi = newRois{r2};
                
                % read warped image
                tmp = roiIndices{r2};
                
                if (gridFacV(1)*gridFacV(2))==1 
                    % if inplane 1x1, then only need to read once!
                    if zIdx==1
                        thisSubROIDir = [sprintf('%s%d_%d_%d',thisDir,tmp(1),tmp(2),tmp(3)) filesep];
                        fname = [thisSubROIDir 'result.mha'];
                        warpMHA = mha_read_volume(fname);
                    end
                else
                    % in-plane subdivision, need to read for each z
                    thisSubROIDir = [sprintf('%s%d_%d_%d',thisDir,tmp(1),tmp(2),tmp(3)) filesep];
                    fname = [thisSubROIDir 'result.mha'];
                    warpMHA = mha_read_volume(fname);
                end
                if  r2==1 && zIdx==1
                    disp(['itksnap -o ' regiBaseDir 'volumes/result.mha -g ' fname ' &'])
                end
                
                % determine HR pixel coordinates
                minIdx = points_to_index( minCoord(:,r2),(ps_uct_1x*res/ps_fac));
                maxIdx = points_to_index( maxCoord(:,r2),(ps_uct_1x*res/ps_fac));
                if r1==1
                    % save starting index, to keep canvas to minimum
                    startMinIdx = minIdx;
                    % initialize canvas etc.
                    %canvasI = zeros(outSize(1:2).*roiIndices{end}(1:2));
                    clear canvasI
                end
                if doTesting && r1==1 && r2==1 && zIdx==1
                    disp(['x ' num2str(r2) ' ' num2str(thisroi(1:2),'%.1f ') 'px, ' ...
                        num2str([minCoord(1,r2) maxCoord(1,r2)],'%.3f ') 'um HR, '  ...
                        num2str([minIdx(1) maxIdx(1)],'%.1f ') 'px HR, res '  ...
                        num2str(outRes(1)) 'um']) %, size ' num2str(outSize(1)) 'px'])
                end
                % full slice
                canvasI([minIdx(1):maxIdx(1)]-startMinIdx(1)+1,[minIdx(2):maxIdx(2)]-startMinIdx(2)+1)=warpMHA(:,:,zIdx);
            end
            
            % output as Tiff slice
            [csy,csx] = size(canvasI);
            tagstruc.ImageLength = csy; % y
            tagstruc.ImageWidth = csx; % x
            tagstruc.BitsPerSample = 16;
            tagstruc.SamplesPerPixel = 1;
            tagstruc.Compression = Tiff.Compression.None;
            tagstruc.SampleFormat = Tiff.SampleFormat.Int;
            tagstruc.Photometric = Tiff.Photometric.MinIsBlack;
            tagstruc.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            
            t = Tiff([outDir 'reco_b' num2str(res) '_warped_' num2str(zV(zIdx),'%04d') '.tif'], 'w');
            t.setTag(tagstruc); t.write(int16(canvasI)); t.close();
            
            if doTesting && zIdx==zTest && gridZ==1 && r==length(roiNameList)
                % show specific canvas for comparison
                
                % stitched result
                figure(8)
                imagesc(canvasI)
                axis equal tight
                title(['stitched HR, res. ' num2str(res*ps_uct_1x) 'um'])
                colormap('gray')
                
                % compare to large roi in initial step, i.e. not stitched
                figure(5)
                imagesc(warpSmallRoiMHA(:,:,zTest))
                title(['not stitched HR, res. ' num2str(res*ps_uct_1x) 'um'])
                axis equal tight
                colormap('gray')
                
                figure(51)
                diffI=double(warpSmallRoiMHA(:,:,zTest))-double(canvasI);
                fidx=find(abs(diffI(:))>0);
                noDiff=length(fidx);
                imagesc(diffI) %,[-0.1 0.1])
                if noDiff>0 && noDiff<100
                    [i1,i2] = ind2sub(size(warpSmallRoiMHA(:,:,zTest)),fidx);
                    hold on
                    % mark positions which are not zero
                    imcontour(abs(diffI)>0,1,'r')
                    plot(i1,i2,'gx','MarkerSize',5)
                    hold off
                end
                axis equal tight
                colormap('gray')
                meanAbsDiff=mean(mean(abs(diffI)));
                title(['diff not stitched - stitched HR, res. ' num2str(res*ps_uct_1x) 'um, meanAbsDiff ' num2str(meanAbsDiff,'%.4f') ' #' num2str(noDiff)])
                colorbar
                
                % downsample to size of fixed image resolution
                figure(6)
                tmpI=imresize3(warpSmallRoiMHA,res*ps_uct_1x/ps_atlas_regi);
                imagesc(tmpI(:,:,zTest))
                title(['downsample not stitched HR, res. ' num2str(ps_atlas_regi) 'um'])
                axis equal tight
                colorbar
                colormap('gray')
                
                % TESTING, show roi of warped image
                tmpLR=double(warped_full_lr(roi(1):roi(2)-1,roi(3):roi(4)-1,roi(5):roi(6)-1));
                figure(7)
                imagesc(tmpLR(:,:,zTest))
                title(['LR warped image, res. ' num2str(ps_atlas_regi) 'um'])
                axis equal tight
                colorbar
                colormap('gray')
                
                % show middle patch of this z
                gStr = num2str(floor(gridFacV(3)/2));
                disp(['itksnap -o ' regiBaseDir 'volumes/result.mha -g ' thisDir gStr '_' gStr '_' num2str(tmp(3)) '/result.mha'])
                %keyboard
            end
            
        end
    end
end
% output runtimes
part2Time = toc(startTiffTime);
disp('Runtime')
disp(['warping subregions ' num2str(part1Time,'%.2f') ' sec, ' num2str(part1Time/60,'%.2f') ' min, ' num2str(part1Time/3600,'%.2f') ' h'])
disp(['writing tiffs      ' num2str(part2Time,'%.2f') ' sec, ' num2str(part2Time/60,'%.2f') ' min, ' num2str(part2Time/3600,'%.2f') ' h'])
disp(['total              ' num2str(part1Time+part2Time,'%.2f') ' sec, ' num2str((part1Time+part2Time)/60,'%.2f') ' min, ' num2str((part1Time+part2Time)/3600,'%.2f') ' h'])
diary off
