function distForEntLvl(lvl)

p1v = 0:0.001:1;
p2v = zeros(size(p1v));
p3v = zeros(size(p1v));
valid = false(size(p1v));
for ii=1:1:length(p1v)
    p1 = p1v(ii);
    try
        p2v(ii) = fzero(@(x)(-p1*zlog2(p1)-x*zlog2(x)-(1-p1-x)*zlog2(1-p1-x) - lvl), 0.5);
        p3v(ii) = 1 - p1 - p2v(ii);
    catch ME
        ME.identifier
        valid(ii) = false;
    end
end

p1v = p1v(valid); p2v = p2v(valid); p3v = p3v(valid);

figure; plot(p1v);
figure; plot(p2v);
figure; plot(p3v);

    function v = zlog2(x)
        if x <= 0
            v = 0;
        else
            v = log2(x);
        end
    end

end