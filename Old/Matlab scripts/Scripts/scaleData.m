function [result,stds]=scaleData(X, direction)
%centerData.m Scales your data matrix.
% Inputs:
%   X   : Raw input matrix of data
%
%   direction: 1, 2 or 3 if the matrix is 3-way. Specifies the mode across
%   which the operation should be performed.
%
% Outputs:
%   result: your scaled matrix of data.
%   stds  : the standard deviations within slabs of the specified direction

if nargin<2
    error('The inputs X and direction must be given.')
end

dim = size(X);
dim_reduced = dim;                      % Copy to make the other dimensions accessible.
dim_reduced(direction) = [];            % Keeps only the dimensions of the modes not being scaled within.
direction_reduced = 1:size(size(X), 2); % For rotation of the array
direction_reduced(direction) = [];

Xperm = permute(X, [direction direction_reduced]);
Xmodified = reshape(Xperm, [dim(direction) prod(dim_reduced)]);

result = zeros(dim(direction), prod(dim_reduced));
stds = zeros(dim(direction), 1);
for i=1:dim(direction)
    dataVector = Xmodified(i,:);
    result(i,:) = dataVector / std(dataVector, 'omitnan');
    stds(i) = std(dataVector, 'omitnan');
end

result = reshape(result, size(Xperm));

new_permutation = zeros(1,size(size(X), 2));
for i=1:size(size(X), 2)
    new_permutation(size(X)==size(result,i)) = i;
end
result = permute(result, new_permutation);