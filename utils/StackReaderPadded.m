function [vol] = StackReaderPadded(readDir,roi)
zstep = 1;

[sx,sy,sz] = StackGetSize(readDir);
loaddir = [basedir sname filesep 'imstack' filesep '*.tif'];
dirinfo = dir(loaddir);

% re-define in terms of pixels defined in dataset
tmp_roi = [max(roi(1),1),min(roi(2),sx),max(roi(3),1),min(roi(4),sy),max(roi(5),1),min(roi(6),sz)];

zlist = tmp_roi(5):zstep:tmp_roi(6);
vol = zeros(tmp_roi(4)-tmp_roi(3)+1,tmp_roi(2)-tmp_roi(1)+1,length(zlist),'uint16');
for i = 1:length(zlist)
    tmpim = imread([dirinfo(zlist(i)).folder filesep dirinfo(zlist(i)).name],...
        'PixelRegion',{[tmp_roi(3),tmp_roi(4)],[tmp_roi(1),tmp_roi(2)]});
    vol(:,:,i) = tmpim;
end

% put zeros for pixels undefined in dataset
if ne(tmp_roi(1),roi(1))
    vol = padarray(vol,[0,tmp_roi(1)-roi(1),0],0,'pre');
end
if ne(tmp_roi(2),roi(2))
    vol = padarray(vol,[0,roi(2)-tmp_roi(2),0],0,'post');  
end

if ne(tmp_roi(3),roi(3))
    vol = padarray(vol,[tmp_roi(3)-roi(3),0,0],0,'pre');
end
if ne(tmp_roi(4),roi(4))
    vol = padarray(vol,[roi(4)-tmp_roi(4),0,0],0,'post');  
end

if ne(tmp_roi(5),roi(5))
    vol = padarray(vol,[0,0,tmp_roi(5)-roi(5)],0,'pre');
end
if ne(tmp_roi(6),roi(6))
    vol = padarray(vol,[0,0,roi(6)-tmp_roi(6)],0,'post');  
end

end
