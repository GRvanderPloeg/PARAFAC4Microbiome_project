% model simulated data
addpath(".\Matlab scripts\Scripts\"); % own scripts
addpath(".\Matlab scripts\N-way toolbox\"); % from Rasmus Bro

% Load data
df_raw = readmatrix("simData.csv", Delimiter=",");
df = df_raw(:,1:J);

subjectMeta = readmatrix("subjectMetadata.csv", Delimiter=",", OutputType="string");
featureMeta = readmatrix("featureMetadata.csv", Delimiter=",", OutputType="string");
subjectLoadings = readmatrix("subjectLoadings.csv", Delimiter=",");
featureLoadings = readmatrix("featureLoadings.csv", Delimiter=",");
timeLoadings = readmatrix("timeLoadings.csv", Delimiter=",");

I = max(str2double(subjectMeta(:,1)));
J = max(str2double(featureMeta(:,1)));
K = size(timeLoadings,2);
numComponents = size(timeLoadings,1);

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
pfac = parafac(df_cnt_scl, numComponents);

% Diagnostic plots for PCA
subplot(2,2,1); scatter(coeff(:,1), coeff(:,2)); title("PCA loadings");
subplot(2,2,2); scatter(score(:,1), score(:,2)); title("PCA scores");
subplot(2,2,3); bar(subjectLoadings); title("Real subject mode");
subplot(2,2,4); bar(featureLoadings); title("Real feature mode");

% Diagnostic plots for PARAFAC
subplot(2,3,1); bar(pfac{1}); title("PARAFAC subject mode");
subplot(2,3,2); bar(pfac{2}); title("PARAFAC feature mode");
subplot(2,3,3); plot(1:9, pfac{3}(:,1)); hold on; plot(1:9, pfac{3}(:,2)); title("PARAFAC time modes");
subplot(2,3,4); bar(subjectLoadings); title("Real subject mode");
subplot(2,3,5); bar(featureLoadings); title("Real feature mode");
subplot(2,3,6); plot(1:9, timeLoadings(1,:)); hold on; plot(1:9, timeLoadings(2,:)); title("Real time modes");