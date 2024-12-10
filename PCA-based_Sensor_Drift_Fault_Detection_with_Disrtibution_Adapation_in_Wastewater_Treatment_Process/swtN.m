% Performs an n-dimensional multilevel discrete wavelet transform of X

% TODO: Returns the coefficients of the 2^n group, approximating the first group (LLL, LLH, LHL, LHH ......) .
% Currently returns approximate coefficients.

function [X] = swtN(X, level, varargin)
    if (ndims(X) == 2)
       
        [X, tmp, tmp, tmp] = swt2(X, level, varargin{:});
        return;
    end

    for dimidx = 1:ndims(X)
     
        xsize = size(X);
        X = unfold(X,2)';

        workingslice = [];
        for modeidx = 1:size(X, 1)  
            [wtemp, tmp] = swt(X(modeidx, :), level, varargin{:});

            if (isempty(workingslice))
                workingslice = zeros(size(X, 1), size(wtemp, 2));
            end

            workingslice(modeidx, :) = wtemp;
        end

       
        xsize(2) = size(workingslice, 2);

        X = shiftdim(fold(workingslice', 2, xsize), 1);
    end
end
