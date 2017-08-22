function directions = oriFieldDirections(ofld)
% Determines the orientation angle of each orientation field point.  This 
% isn't a true vector field; elements pointing in opposite directions
% are equivalent, and so the directions will be in the range [0, pi).
% INPUTS:
%   ofld: MxNx2 orientation field for which the angles are to be determined
% OUTPUT:
%   directions: MxN matrix of angles of each orientation field point.
%       Angles are in radians and are in the range [0, pi).

directions = mod(atan(ofld(:, :, 2) ./ ofld(:, :, 1)), pi);
directions(ofld(:, :, 1) == 0) = 0;

end