function dist = angleDist(aFrom, aTo)
% arc length between aFrom and aTo on a circle (i.e., with a wraparound at
% 2*pi)

wrapArcs = aFrom > aTo;
dist = (aTo - aFrom);
dist(wrapArcs) = dist(wrapArcs) + 2*pi;

% if aFrom < aTo
%     dist = aTo - aFrom;
% else % wraparound
%     dist = aTo - aFrom + 2*pi;
% end

end