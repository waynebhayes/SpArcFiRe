classdef PriorityQueue < handle
    properties (Hidden, SetAccess = private)
        keys = [];
        vals = [];
%         locs = [];
%         origLocs = [];
        nElt = 0;
        locKeys = [];
        origLocKeys = [];
        loc2idx = [];
    end
    methods
        function PQ = PriorityQueue(keys, vals)
            if nargin == 1
                error('if keys specified, values must also be specified');
            elseif nargin == 2
                if length(keys) ~= length(vals)
                    size(keys)
                    size(vals)
                    error('key and value arrays must have the same length');
                end
                if ~isempty(keys)
                    if ~isvector(keys) || ~isvector(vals)
                        error('keys and values must be expressed as a vector');
                    end
                    if size(keys, 2) > 1
                        keys = keys';
                    end
                    if size(vals, 2) > 1
                        vals = vals';
                    end
                    PQ.nElt = numel(keys);
                    PQ.keys = keys;
                    PQ.vals = vals;
%                     PQ.locs = Ref2Num(1:1:PQ.nElt);
%                     PQ.origLocs = PQ.locs;
                    PQ.locKeys = int32(1:1:PQ.nElt);
                    PQ.origLocKeys = PQ.locKeys;
                    PQ.loc2idx = containers.Map(num2cell(PQ.locKeys), num2cell(PQ.locKeys));
                    heapify(PQ);
                end
            end
%             PQ.keys
%             PQ.vals
%             size(PQ.keys)
%             size(PQ.vals)
%             assert(isempty(PQ.keys) || isvector(PQ.keys));
%             assert(size(PQ.keys, 2) <= 1);
%             assert(isempty(PQ.vals) || isvector(PQ.vals));
%             assert(size(PQ.vals, 2) <= 1);
        end
        function pqLocs = getStartLocs(PQ)
            % Locators for the keys and values given at priority queue
            % instantiation.  These locator objects should not be modified
            % externally.
%             pqLocs = PQ.origLocs;
            pqLocs = PQ.origLocKeys;
        end
        function [key, val] = peekMax(PQ)
            if PQ.nElt < 1
                error('attempt to access element of empty priority queue');
            end
            key = PQ.keys(1);
            val = PQ.vals(1);
        end
        function [key, val] = takeMax(PQ)
            [key, val] = peekMax(PQ);
            removeEltByIdx(PQ, 1);
        end
        function changeKeyAtLoc(PQ, newKey, loc)
%             idx = loc.val;
            idx = PQ.loc2idx(loc);
            PQ.keys(idx) = newKey;
            siftUp(PQ, idx);
            siftDown(PQ, idx);
        end
        function removeEltAtLoc(PQ, loc, dbg)
            if nargin > 2 && dbg
                idx = PQ.loc2idx(loc)
            end
            PQ.removeEltByIdx(PQ.loc2idx(loc));
        end
        function epty = isEmpty(PQ)
%             assert(PQ.nElt >= 0);
%             assert((numel(PQ.keys) - sum(isnan(PQ.keys))) == PQ.nElt);
            epty = (PQ.nElt < 1);
        end
        function makeEmpty(PQ)
            PQ.keys = [];
            PQ.vals = [];
            PQ.nElt = 0;
            PQ.locKeys = [];
            PQ.origLocKeys = [];
            PQ.loc2idx = [];
        end
        function v = getVals(PQ)
            v = PQ.vals(1:PQ.nElt);
        end
        function k = dbgGetKeys(PQ)
            k = PQ.keys(1:PQ.nElt);
        end
        function display(PQ)
            fprintf('Priority queue of size %d\n', PQ.nElt);
        end
        function dbgDisplay(PQ)
            fprintf('%d elements\n', PQ.nElt);
            fprintf('keys: %s\n', mat2str(PQ.keys(1:PQ.nElt), 6));
            fprintf('vals: %s\n', mat2str(PQ.vals(1:PQ.nElt), 6));
        end
        function v = dbgGetValAtLoc(PQ, loc)
            v = PQ.vals(PQ.loc2idx(loc));
        end
    end
    methods (Access = private)
        function siftUp(PQ, idx)
            pIdx = floor(idx/2);
            while( pIdx > 0 && (PQ.keys(pIdx) < PQ.keys(idx)) )
                PQ.keys([pIdx idx]) = PQ.keys([idx pIdx]);
                PQ.vals([pIdx idx]) = PQ.vals([idx pIdx]);
%                 PQ.locs([pIdx idx]) = PQ.locs([idx pIdx]);
%                 PQ.locs(pIdx).val = pIdx;
%                 PQ.locs(idx).val = idx;
                PQ.locKeys([pIdx idx]) = PQ.locKeys([idx pIdx]);
                PQ.loc2idx(PQ.locKeys(idx)) = idx;
                PQ.loc2idx(PQ.locKeys(pIdx)) = pIdx;
                idx = pIdx;
                pIdx = floor(idx/2);
            end
        end
        function siftDown(PQ, idx)
%             fprintf('siftDown\n');
            while idx <= PQ.nElt
%                 fprintf('\tidx = %d\n', idx);
                lIdx = 2*idx; rIdx = lIdx + 1;
                if ( rIdx <= PQ.nElt && (PQ.keys(rIdx) > PQ.keys(lIdx)) )
                    lgrIdx = rIdx;
                else
                    lgrIdx = lIdx;
                end
                if ( lgrIdx <= PQ.nElt && (PQ.keys(lgrIdx) > PQ.keys(idx)) )
                    PQ.keys([lgrIdx idx]) = PQ.keys([idx lgrIdx]);
                    PQ.vals([lgrIdx idx]) = PQ.vals([idx lgrIdx]);
%                     PQ.locs([lgrIdx idx]) = PQ.locs([idx lgrIdx]);
%                     PQ.locs(lgrIdx).val = lgrIdx;
%                     PQ.locs(idx).val = idx;
                    PQ.locKeys([lgrIdx idx]) = PQ.locKeys([idx lgrIdx]);
                    PQ.loc2idx(PQ.locKeys(idx)) = idx;
                    PQ.loc2idx(PQ.locKeys(lgrIdx)) = lgrIdx;
                    idx = lgrIdx;
                else
                    idx = inf;
                end
            end
        end
        function heapify(PQ)
            for idx = floor(PQ.nElt / 2) : -1 : 1
                siftDown(PQ, idx);
            end
        end
        function removeEltByIdx(PQ, idx)
            lastIdx = PQ.nElt;
            PQ.keys(idx) = PQ.keys(lastIdx);
            PQ.vals(idx) = PQ.vals(lastIdx);
%             PQ.locs(idx) = PQ.locs(lastIdx);
            PQ.locKeys(idx) = PQ.locKeys(lastIdx);
            PQ.loc2idx(PQ.locKeys(idx)) = idx;
            PQ.keys(lastIdx) = NaN;
            PQ.vals(lastIdx) = NaN;
%             PQ.locs(lastIdx).val = NaN;
%             PQ.locs(lastIdx) = [];
            PQ.locKeys(lastIdx) = -1;
            PQ.nElt = PQ.nElt - 1;
            siftUp(PQ, idx);
            siftDown(PQ, idx);
        end
    end
    methods (Static)
        function [PQm, size, locs] = merge(PQ1, PQ2, keep1, keep2)
            if nargin == 3
                error('if keep1 specified, keep2 must also be specified');
            end
            keys1 = PQ1.keys(1:PQ1.nElt);
            vals1 = PQ1.vals(1:PQ1.nElt);
            keys2 = PQ2.keys(1:PQ2.nElt);
            vals2 = PQ2.vals(1:PQ2.nElt);
%             if ~isempty(rmVal1)
%                 keep1 = (vals1 ~= rmVal1);
%                 keys1 = keys1(keep1); vals1 = vals1(keep1);
%                 keep2 = (vals2 ~= rmVal2);
%                 keys2 = keys2(keep2); vals2 = vals2(keep2);
%             end
            if ~isempty(keep1)
                keys1 = keys1(keep1); vals1 = vals1(keep1);
            end
            if ~isempty(keep2)
                keys2 = keys2(keep2); vals2 = vals2(keep2);
            end
            PQm = PriorityQueue([keys1; keys2], [vals1; vals2]);
            size = PQm.nElt;
            if nargout > 2
                locs = PQm.getStartLocs();
            end
        end
    end
end