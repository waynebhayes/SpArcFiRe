function mtx = polXyMtx(sz, ctrX, ctrY, polXyFxn, params, thBounds)

[x, y] = polXyFxn(params, thBounds);

xRange = [1 sz(2)] - ctrX;
yRange = [1 sz(1)] - ctrY;
inRange = (x >= min(xRange)) & (x < max(xRange)) & (y >= min(yRange)) & (y < max(yRange));
x = x(inRange);
y = y(inRange);

x = round(x + ctrX);
y =  round(y + ctrY);
mtx = zeros(sz);
mtx(sub2ind(sz, y, x)) = 1;
mtx = flipud(mtx);

end

