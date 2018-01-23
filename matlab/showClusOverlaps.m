function showClusOverlaps(clusMtxList1, clusMtxList2)

imgSize = size(clusMtxList1(:, :, 1));

cImg1 = showClustersFromMtxs(clusMtxList1, imgSize);
cImg2 = showClustersFromMtxs(clusMtxList2, imgSize);
figure;
subplot(1, 2, 1); imshow(cImg1); title('waveband 1');
subplot(1, 2, 2); imshow(cImg2); title('waveband 2');

cbndImg1 = displayClusterOverlay(zeros(imgSize), clusMtxList1);
cbndImg2 = displayClusterOverlay(zeros(imgSize), clusMtxList2);

figure; imshow(cbndImg2(:, :, [2 1 3]) + cbndImg1);

end