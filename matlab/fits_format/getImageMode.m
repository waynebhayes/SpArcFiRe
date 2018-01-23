function histMode = getImageMode(img, nBins)

imgVec = img(:);

edges = linspace(min(imgVec), max(imgVec), nBins);
[n, binIdxs] = histc(imgVec, edges);

% [n, xout] = hist(img(:), nBins);

% figure; bar(xout, n);

[mV, mI] = max(n);

% histMode = xout(mI);

% histMode = mean(edges(mI:mI+1))
histMode = mean(imgVec(binIdxs == mI));

figure; hist(imgVec(binIdxs == mI), 100);

end