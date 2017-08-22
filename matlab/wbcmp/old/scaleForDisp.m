
function scImg = scaleForDisp(img)
    lims = quantile(nonzeros(img), [0.1, 0.9]);
    scImg = img;
    scImg = scImg - lims(1);
    scImg = scImg ./ (lims(2) - lims(1));
    scImg(scImg < 0) = 0;
    scImg(scImg > 1) = 1;
end