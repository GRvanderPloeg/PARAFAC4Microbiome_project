% model simulated data
addpath(".\Matlab scripts\Scripts\"); % own scripts
addpath(".\Matlab scripts\N-way toolbox\"); % from Rasmus Bro

I = 20; % number of subjects
J = 100; % number of features
K = 9; % number of timepoints

df_raw = readmatrix("simData.csv", Delimiter=",");
df = df_raw(:,1:J);

% processing
df_clr = transformCLR(df);
df_cube = reshape(df_clr, I, J, K);
[df_cnt, df_means] = centerData(df_cube, 1);
[df_cnt_scl, df_stds] = scaleData(df_cnt, 2);

df_cnt_scl_flat = df_clr;
for i=1:size(df_clr,2)
    df_cnt_scl_flat(:,i) = (df_clr(:,i) - mean(df_clr(:,i))) /  std(df_clr(:,i));
end

% modelling
[coeff, score, latent, tsquared, explained] = pca(df_cnt_scl_flat, 'Centered',false);
pfac = parafac(df_cnt_scl, 2);

subplot(2,3,1); scatter(coeff(:,1), coeff(:,2));
subplot(2,3,2); scatter(score(:,1), score(:,2));
subplot(2,3,3); scatter(pfac{1}(:,1), pfac{1}(:,1));
subplot(2,3,4); scatter(pfac{2}(:,1), pfac{2}(:,1));
subplot(2,3,5); plot(1:9, pfac{3}(:,1));
subplot(2,3,6); plot(1:9, pfac{3}(:,2));