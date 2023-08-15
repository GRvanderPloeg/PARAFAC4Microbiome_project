function [keepIndex,sparsityPerFeature,variationPerFeature]=removeColumns(Xraw, sparsityThreshold, variationThreshold)
[I,J,K] = size(Xraw);

dummy = reshape(Xraw, I, J*K)';
Xlong = reshape(dummy, J, I*K)';

totalSSQ = sumsqr(Xlong);
sparsityPerFeature = zeros(J,1);
variationPerFeature = zeros(J,1);
keepIndex = false(1,J); % Logical vector that keeps track of which features should be kept.

for j=1:J
    dataVector = Xlong(:,j);
    sparsity = sum(dataVector==0)/(I*K) * 100;
    variation = sumsqr(dataVector) / totalSSQ * 100;
    sparsityPerFeature(j) = sparsity;
    variationPerFeature(j) = variation;

    if(sparsity <= sparsityThreshold && variation >= variationThreshold)
        keepIndex(j) = true;
    end
end

keepIndex = keepIndex';