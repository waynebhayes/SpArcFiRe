classdef Ref2Num < handle
% old experiment in using OOP/references/priority queues in Matlab. Matlab
% is slow at this; do not use.
    properties
        val = NaN;
    end
    methods
        function R2N = Ref2Num(arr)
            if nargin > 0
                if isscalar(arr)
                    R2N.val = arr;
                else
                    n = length(arr);
                    R2N(n) = Ref2Num();
                    for ii=1:1:n
                        R2N(ii) = Ref2Num(arr(ii));
                    end   
                end
            end
        end
        function display(R2N)
            fprintf('ref to num with val = %f\n', R2N.val);
        end
        function isEq = eq(a, b)
            if ( (length(a) ~= length(b)) && ~isscalar(a) && ~isscalar(b) )
                error('nonscalar a and b must have the same length');
            end
            if isscalar(a) && isscalar(b)
                isEq = (a.val == b.val);
                return
            end
            % avoid copying arrays - make more concise if Matlab does this
            % automatically
            if ~isscalar(a) && isscalar(b)
                isEq = false(size(a));
                for ii=1:1:length(a)
                    isEq(ii) = (a(ii) == b);
                end
            elseif isscalar(a) && ~isscalar(b)
                isEq = false(size(b));
                for ii=1:1:length(b)
                    isEq(ii) = (a == b(ii));
                end
            else
                isEq = false(size(a));
                for ii=1:1:length(a)
                    isEq(ii) = (a(ii) == b(ii));
                end
            end
%             isEq = ( arrayfun(@(x)(x.val), a) == arrayfun(@(x)(x.val), b) );
        end
        function isNe = ne(a, b)
            isNe = ~eq(a, b);
        end
    end
%     methods (Static)
%         function c = copy(R2N)
%             c = R2N;
%             if ~isscalar(R2N)
%                 for ii=1:1:numel(R2N)
%                     c(ii) = R2N(ii);
%                 end
%             end
%         end
%     end
end
        