function [df, modeI, modeJ, modeK]=dataFromVectorsPlusNoise(numComponents, I, J, K, Istd, Jstd, Kstd, Imean, Jmean, Kmean, Inoise_std, Jnoise_std, Knoise_std, overallNoise_std)
[df, modeI, modeJ, modeK] = dataFromVectorsNoNoise(numComponents, I, J, K, Istd, Jstd, Kstd, Imean, Jmean, Kmean);

% Generate and add noise
id_noise = generateNormalRandomValues(I, Inoise_std, 0);
feature_noise = generateNormalRandomValues(J, Jnoise_std, 0);
time_noise = generateNormalRandomValues(K, Knoise_std, 0);
overall_noise = overallNoise_std.*randn(I,J,K);

for i=1:I
    df(i,:,:) = df(i,:,:) + id_noise(i);
end
for j=1:J
    df(:,j,:) = df(:,j,:) + feature_noise(j);
end
for k=1:K
    df(:,:,k) = df(:,:,k) + time_noise(k);
end
df = df + overall_noise;