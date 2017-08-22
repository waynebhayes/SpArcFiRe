function plot_pang_cmp(gxyParams, cmpPangs, qmeasIdx, qmeasName)

if nargin < 4 || isempty(qmeasName)
    qmeasName = '';
end

arc_params = gxyParams.arc_params;

qmeas = arc_params(:, qmeasIdx);
measPangs = arc_params(:, 2);
figure; scatter(abs(measPangs), qmeas);
xlabel('pitch angle (degrees)');
ylabel(qmeasName);

for ii=1:1:length(cmpPangs)
    vline(cmpPangs(ii));
end
vline(mean(abs(measPangs)), 'b:')
vline(sum(abs(measPangs) .* qmeas) ./ sum(qmeas), 'b-')

end