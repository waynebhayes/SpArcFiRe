function y = cmhf(x)

y = (2 / sqrt(3)) * (pi^(-1/4)) * (sqrt(pi) * (1 - x.^2) .* exp((-1/2)*x.^2) - ...
    (sqrt(2)*1i*x + sqrt(pi)*erfz((1i/sqrt(2)) .* x) .* (1 - x.^2) .* exp((-1/2)*x.^2)));

end