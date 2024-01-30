function [df, modeI, modeJ, modeK]=dataFromVectorsNoNoise(numComponents, I, J, K, Istd, Jstd, Kstd, Imean, Jmean, Kmean)
modeI = zeros(I, numComponents);
modeJ = zeros(J, numComponents);
modeK = zeros(K, numComponents);

for n=1:numComponents
    Ivector = generateNormalRandomValues(I, Istd, Imean);
    Ivector = Ivector ./ sqrt(sum(Ivector.^2)); % normalize to length 1.
    modeI(:,n) = Ivector;

    Jvector = generateNormalRandomValues(J, Jstd, Jmean);
    Jvector = Jvector ./ sqrt(sum(Jvector.^2)); % normalize to length 1.
    modeJ(:,n) = Jvector;

    Kvector = generateNormalRandomValues(K, Kstd, Kmean);
    Kvector = Kvector ./ sqrt(sum(Kvector.^2)); % normalize to length 1.
    modeK(:,n) = Kvector;
end

df = modeI*krb(modeK,modeJ)';
df = reshape(df, I, J, K);