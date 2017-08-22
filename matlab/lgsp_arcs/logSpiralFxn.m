function rho = logSpiralFxn(theta, params)
thOff = params(1); pitchAngle = params(2); initRadius = params(3);

rho = initRadius * exp(-pitchAngle * mod((theta - thOff), 2*pi));

end