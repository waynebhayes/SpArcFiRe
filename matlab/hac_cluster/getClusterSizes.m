function sizes = getClusterSizes(clusters)

sizes = arrayfun(@(x)(length(x.mergedPts)), clusters);

end