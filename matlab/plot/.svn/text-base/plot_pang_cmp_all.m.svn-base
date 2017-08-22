function plot_pang_cmp_all(gxyParams, cmpPangs)

gxyName = gxyParams.name;

plot_pang_cmp(gxyParams, cmpPangs, 6, 'arc lengths');
saveas(gcf, [gxyName '_arcLengths.png']);

plot_pang_cmp(gxyParams, cmpPangs, 7, 'pixel counts');
saveas(gcf, [gxyName '_pixCounts.png']);

plot_pang_cmp(gxyParams, cmpPangs, 10, 'average brightness');
saveas(gcf, [gxyName '_avgBrt.png']);

plot_pang_cmp(gxyParams, cmpPangs, 11, 'clus-scores');
saveas(gcf, [gxyName '_clusScores.png']);

plot_pang_cmp(gxyParams, cmpPangs, 12, 'usm-clus-scores');
saveas(gcf, [gxyName '_clusScoresUsm.png']);

end