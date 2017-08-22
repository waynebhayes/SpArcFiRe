function phase = getPhase(img, angle)
% old experiment, do not use

filtLight = genOriFilterFxn(angle, [], 1);
filtDark = genOriFilterFxn(angle, [], -1);

fimgL = conv2(img, filtLight, 'same');
fimgD = conv2(img, filtDark, 'same');

stgths = sqrt(fimgL.^2 + fimgD.^2);
figure; imagesc(fimgL); axis image
figure; imagesc(fimgD); axis image
figure; imagesc(stgths); axis image

phase = atan2(fimgL, abs(fimgD));

end