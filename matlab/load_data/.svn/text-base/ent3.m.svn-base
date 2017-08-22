function ent = ent3(p1, p2, p3)
% entropy for a probability density function with three possible outcomes,
% with probabilities p1, p2, and p3

if nargin < 3
    p3 = 1 - (p1 + p2);
end

pvs = [p1(:) p2(:) p3(:)];

ent = pvs .* log2(pvs);
ent(pvs == 0) = 0;
ent = -sum(ent, 2);

ent = reshape(ent, size(p1));

end