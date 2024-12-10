% Matrixes (unfolds) a tensor along a specified pattern.
function [mat] = unfold(tensor, mode)
    othermodes = setdiff(1:ndims(tensor), mode);
    mat = reshape(permute(tensor, [mode othermodes]), [size(tensor, mode) (numel(tensor) ./ size(tensor, mode))]);
end
