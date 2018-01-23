function dofld = deNoiseOrientationField(ofld)
% Performs the orientation field de-noising described in the PhD thesis
% "Inferring Galaxy Morphology Through Texture Analysis" (K. Au 2006).
% INPUTS:
%   ofld: orientation field without de-noising (but already merged from 3
%       resolution levels if the process in the thesis is followed)
% OUTPUTS:
%   dofld: de-noised orientation field

if length(size(ofld)) ~= 3 || size(ofld, 3) ~= 2
    error('input orientation field should be a m x n x 2 matrix')
end

% We do this to keep the orientation field vectors from all being zeroed
% out during de-noising.  This is the only known deviation from the process
% described in the thesis.
% UPDATE: not needed anymore, after we fixed the sampling points for the
% orientation filters
% mean(mean(oriFieldStrengths(ofld)))
% ofld = ofld / mean(mean(oriFieldStrengths(ofld)));
% ofld = ofld * 10;

% The neighbor distance is 5 in the thesis; if needed, this could be
% changed here or made a parameter
nbrDist = 5;
dofld = zeros(size(ofld));

% TODO: vectorize this without losing generality
oriRows = size(ofld, 1);
oriCols = size(ofld, 2);
% warning('using different orientation-field de-noising');
for r=1:1:oriRows
    for c=1:1:oriCols
        curVec = squeeze(ofld(r, c, :));
        normCurVec = norm(curVec);
        nbrVecs = [];
        if r - nbrDist >= 1 && c - nbrDist >= 1
            nbrVec = ofld(r - nbrDist, c - nbrDist, :);
            nbrVecs(end+1, 1:2) = nbrVec;
        end
        if r + nbrDist <= oriRows && c - nbrDist >= 1
            nbrVec = ofld(r + nbrDist, c - nbrDist, :);
            nbrVecs(end+1, 1:2) = nbrVec;
        end
        if r - nbrDist >= 1 && c + nbrDist <= oriCols
            nbrVec = ofld(r - nbrDist, c + nbrDist, :);
            nbrVecs(end+1, 1:2) = nbrVec;
        end
        if r + nbrDist <= oriRows && c + nbrDist <= oriCols
            nbrVec = ofld(r + nbrDist, c + nbrDist, :);
            nbrVecs(end+1, 1:2) = nbrVec;
        end
        nbrSims = zeros(size(nbrVecs, 1), 1);
        
        subtrAmt = cos(pi/4);
        for nbr=1:1:size(nbrVecs, 1)
            nbrSims(nbr) = ...
                max(abs(nbrVecs(nbr, 1:2) * curVec) - subtrAmt, 0) / ...
                (normCurVec * norm(nbrVecs(nbr, 1:2)));
            
%             nbrSims(nbr) = ...
%                 max(abs(nbrVecs(nbr, 1:2) * curVec) / (normCurVec * norm(nbrVecs(nbr, 1:2))) - subtrAmt, 0) ;
            
%             nbrSims(nbr) = max( (abs(nbrVecs(nbr, 1:2) * curVec) / ...
%                 (normCurVec * norm(nbrVecs(nbr, 1:2)))) - cos(pi/4), 0);
            
%             nbrSims(nbr) = ...
%                 max(abs(nbrVecs(nbr, 1:2) * curVec) - (normCurVec * norm(nbrVecs(nbr, 1:2)) * subtrAmt), 0) / ...
%                 (normCurVec * norm(nbrVecs(nbr, 1:2)));
%             nbrSims(nbr) = ...
%                 abs(nbrVecs(nbr, 1:2) * curVec) / ...
%                 (normCurVec * norm(nbrVecs(nbr, 1:2)));
%             nbrSims(nbr) = nbrSims(nbr) * (acos(nbrSims(nbr)) < subtrAmt || pi - acos(nbrSims(nbr)) < subtrAmt);
        end
        dofld(r, c, :) = (curVec / normCurVec) * median(nbrSims);
%         dofld(r, c, :) = curVec * median(nbrSims);
%         dofld(r, c, :) = curVec;
    end
end
end