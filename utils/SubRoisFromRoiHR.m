function [newRois,newIndices] = SubRoisFromRoiHR(roi,gridFacV)
% create non-overlapping sub-rois from an initial roi
% [newRois,newIndices] = SubRoisFromRoiHR(roi,gridFacV)
%
% it will make a grid [gridFacV(1),gridFacV(2),gridFacV(3)] of rois, each of size
% approximately the size of roi ./ gridFacV
%
% roi is given already in high-resolution pixel indices
% hence split should keep integer values of pixel indices

tmp = diff(roi);
roiSize = tmp([1,3,5]);
tmp = ceil(roiSize./gridFacV);  % ensure integer values

newRois = {};
newIndices = {};
c = 1;
for i3 = 1:gridFacV(3)
    for i2 = 1:gridFacV(2)
        for i1 = 1:gridFacV(1)
            ii1 = roi(1)+(i1-1)*tmp(1); % initial index 1st coordinate
            fi1 = roi(1)+ min(i1*tmp(1),roiSize(1)); % final index 1st coordinate
            ii2 = roi(3)+(i2-1)*tmp(2);
            fi2 = roi(3)+min(i2*tmp(2),roiSize(2));
            ii3 = roi(5)+(i3-1)*tmp(3);
            fi3 = roi(5)+ min(i3*tmp(3),roiSize(3));
            thisRoi = [ii1,fi1,ii2,fi2,ii3,fi3];
            newRois{c} = thisRoi;
            newIndices{c} = [i1,i2,i3];
            c = c+1;
        end
    end
end
end

