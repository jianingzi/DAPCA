% % assigns elements to cells and calculates cell counts.
% is intended to be a helper function; don't call it directly unless you have a reason to.
function [counts] = assign_cells(datacells, weights, num_cells)
counts = zeros(num_cells);
idxcell = num2cell(datacells, 1);
countidx = sub2ind(size(counts), idxcell{:});
for countpos = 1:length(countidx)
    counts(countidx(countpos)) = counts(countidx(countpos)) + weights(countpos);
end
% Low order: use sparse. High order: use sptensor.
% counts = sparse(counts);
% Most of the counts will be 0, which saves memory.
end