function df=makeNoise(df,Inoise_std,Jnoise_std,Knoise_std,overallNoise_std)
I = size(df, 1);
J = size(df, 2);
K = size(df, 3);

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