function imgMtxs = loadImagesFromPaths(paths)

validateattributes(paths, {'cell'}, {});

imgMtxs = cell(1, length(paths));
for ii=1:1:length(paths)
    curImg = im2double(imread(paths{ii}));
    imgSz = size(curImg);
    if (length(imgSz) == 3 && imgSz(3) == 3) % color image, make gray
       curImg = rgb2gray(curImg);
    end
    imgMtxs{ii} = curImg;
end

end