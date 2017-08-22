function err = wbCmpLgspErr(lgspFxn, lgspParams, theta, rho, brt)

error(nargchk(4, 5, nargin));

if nargin < 5 || isempty(brt)
    brt = ones(size(theta));
end

err = brt .* sqrt(abs(rho - lgspFxn(theta, lgspParams)));
% err = exp(brt/1000) .* log(abs(rho - lgspFxn(theta, lgspParams))+1);



end