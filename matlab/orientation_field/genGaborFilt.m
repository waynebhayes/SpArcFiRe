function gbf = genGaborFilt(size, wl, angle, sigma)
% experimental

[xc, yc] = meshgrid(linspace(-pi, pi, size), linspace(-pi, pi, size));

x = xc * cos(angle) + yc * sin(angle);
y = -xc * sin(angle) + yc * cos(angle);

figure; imagesc(exp(-(x.^2 + y.^2) / (2*sigma^2))); axis image
figure; imagesc( cos((2*pi*x)/wl) ); axis image

gbf = exp(-(x.^2 + y.^2) / (2*sigma^2)) .* cos((2*pi*x)/wl);

figure; imagesc(gbf); axis image

end