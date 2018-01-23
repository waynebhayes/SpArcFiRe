function err = lgspErr(lgspFxn, lgspParams, theta, rho, brt)
% Fit error for a logarithmic spiral, in the form needed by lsqnonlin

error(nargchk(4, 5, nargin));

if nargin < 5 || isempty(brt)
    brt = ones(size(theta));
end

err = sqrt(brt) .* (rho - lgspFxn(theta, lgspParams));
% err = sqrt(brt) .* sqrt(abs(rho - lgspFxn(theta, lgspParams)));  % TEMP

% figure; hold on; polar(theta, rho, 'g-'); polar(theta, lgspFxn(theta, lgspParams), 'k-'); title(sprintf('%s', mat2str([min(theta) max(theta)]))); hold off
% figure; scatter(theta, err); xlabel('theta'); ylabel('err');
% figure; polar(theta, err, '.'); title('theta-err');

end