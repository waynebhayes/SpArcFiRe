function clusMtxs = wtClustMtxsByFitsBrt(clusMtxs, fitsBrt)

blackLevelQuantile = 0.00;
whiteLevelQuantile = 0.99;

for ii=1:1:size(clusMtxs, 3)
    wts = fitsBrt .* (clusMtxs(:, :, ii) > 0);
    bqL = quantile(wts(wts > 0), blackLevelQuantile);
    wqL = quantile(wts(wts > 0), whiteLevelQuantile);
    wts = wts - bqL;
    wts = wts / wqL;
    wts(wts < 0) = 0;
    wts(wts > 1) = 1;
    clusMtxs(:, :, ii) = wts;
end

figure; imshow(sum(clusMtxs, 3));
figure; imshow(sum(clusMtxs, 3) > 0);

end