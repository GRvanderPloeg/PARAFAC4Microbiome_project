function Xclr=transformCLR(X)
% Centered-log ratio transformation of compositional data
% Zero-handling implemented as suggested in Quinn et al 2019: 10.1093/gigascience/giz107
% Currently limited to two-way datasets
%
% Input:
% X     : input dataset
%
% Output:
% Xclr  : CLR transformed data

[I,J] = size(X);

% Calculate pseudocount as smallest nonzero element and add it to all
vectorizedX = reshape(X, I*J, 1);
pseudocount = min(vectorizedX(vectorizedX>0));
X = X + pseudocount;

Xlong = X;
Xclr = Xlong;

for i=1:I
    dataVector = Xlong(i,:);
    Xclr(i,:) = log(dataVector/geomean(dataVector));
end

