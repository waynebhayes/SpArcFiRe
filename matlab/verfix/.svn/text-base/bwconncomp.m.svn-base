function cc = bwconncomp(bw)
% for MATLAB versions where bwconncomp isn't available

bw = bw > 0;

cc.ImageSize = size(bw);

[lbls, nCC] = bwlabel(bw);
cc.NumObjects = nCC;

idxLists = cell(1, nCC);
for ii=1:1:nCC
    idxLists{ii} = find(lbls == ii);
end

cc.PixelIdxList = idxLists;

end