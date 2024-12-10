% Converts the dataset (observations in rows, features in columns) into a multidimensional grid format.
% Returns a count of each grid cell (and therefore also performs quantisation).
% A scalar can be specified as the same number of cells per dimension, or a vector can be specified as a different number of cells per dimension.

function [counts, datacells] = data2grid(data, weights, num_cells)
   % If the data is not weighted, a weight of 1 is applied to each grid cell.
   if (~exist('weights', 'var') || isempty(weights))
       weights = ones(size(data, 1), 1);
   end

   % per cent to find the extent of the grid.
   featuremin = min(data, [], 1);
   featuremax = max(data, [], 1);

  
    if (isscalar(num_cells))
        num_cells = repmat(num_cells, [1 size(data,2)]);
    end
    
    cellsize = (featuremax + 1 - featuremin) ./ num_cells;
    

    tilesize = size(data);
    tilesize(2) = 1;

    %  Quantise the data to its cell index.
    %  Quantisation doesn't work very well for vectors of quantised values.
    datacells = floor((data - repmat(featuremin, tilesize)) ./ repmat(cellsize, tilesize)) + 1;

    counts = assign_cells(datacells, weights, num_cells);
%    counts = counts ./ sqrt(prod(cellsize)); 

end
