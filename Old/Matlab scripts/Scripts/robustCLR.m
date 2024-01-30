function Xclr=robustCLR(X)
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
Xlong = X;
Xclr = Xlong;

for i=1:I
    dataVector = Xlong(i,:);
    g = geomean(dataVector(dataVector > 0), "omitnan");
    result = log(dataVector/g);
    result(result==-Inf) = NaN;
    Xclr(i,:) = result;
end

