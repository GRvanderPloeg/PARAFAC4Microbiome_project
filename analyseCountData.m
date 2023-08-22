% Analyse count data
% Initialize
addpath(".\Matlab scripts\Scripts\"); % own scripts
addpath(".\Matlab scripts\N-way toolbox\"); % from Rasmus Bro

%%
% Load data
df = readmatrix("tempData.csv", Delimiter=",");

subjectMeta = readmatrix("subjectMetadata.csv", Delimiter=",", OutputType="string");
featureMeta = readmatrix("featureMetadata.csv", Delimiter=",", OutputType="string");
subjectLoadings = readmatrix("subjectLoadings.csv", Delimiter=",");
featureLoadings = readmatrix("featureLoadings.csv", Delimiter=",");
timeLoadings = readmatrix("timeLoadings.csv", Delimiter=",");

I = size(subjectLoadings, 1);
J = size(featureLoadings, 1);
K = size(timeLoadings,1);
numComponents = size(timeLoadings,2);

%%
% Processing
[df_clr, ~] = transformCLR(df);

df_cube = reshape(df_clr, I, K, J);
df_cube = permute(df_cube, [1 3 2]);

[df_cnt, df_means] = centerData(df_cube, 1);
[df_cnt_scl, df_stds] = scaleData(df_cnt, 2);

df_cnt_scl_flat = df_clr;
for i=1:size(df_clr,2)
    df_cnt_scl_flat(:,i) = (df_clr(:,i) - mean(df_clr(:,i))) /  std(df_clr(:,i));
end

%%
% PCA Modelling
[coeff, score, latent, tsquared, explained] = pca(df_cnt_scl_flat, 'Centered',false, "NumComponents",numComponents);
explained(1:5)

%% 
% Initialize options
Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 2;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)

%%
% Do a reporter Parafac + plot it to check what the best model is.
path_start = "./test_run/Figures/";
maxComponents=3;
%days = [-14 0 2 5 9 14 21];
days = 1:K;
numReps=25;
maxIterations=20;

[tongueModels, tongueCons, tongueVarExps, tongueBoots, tongueBootVarExps, tongueTuckers] = quickReport(df_cnt_scl, maxComponents, numReps, maxIterations, subjectMeta, featureMeta, days, 2, "Tongue bootstrapped", path_start+"tongue");

%%
% Dump data so far for later inspection
path_start = "./test_run/Dump/";
dump(tongueModels, tongueCons, tongueVarExps, tongueBoots, tongueBootVarExps, tongueTuckers, path_start, "tongue");

%%
% Pick the number of components that seems appropriate based on the plots
% and pick the model that you like the best.
% Best practices: pick between best Corcondia or best variance explained.
% Alternative: choose a good tucker congruence model.

numFactors_tongue = numComponents;

tongue_choice = find(tongueVarExps{numFactors_tongue}==max(tongueVarExps{numFactors_tongue}));
%tongue_choice = find(tongueCons{numFactors_tongue}==max(tongueCons{numFactors_tongue}));

Tongue_model = pickModel(tongueModels{1,numFactors_tongue}, tongueModels{2,numFactors_tongue}, tongueModels{3,numFactors_tongue}, tongue_choice);

%%
% Save models
model_path = "./test_run/PARAFAC models/";
savePARAFAC(df_cnt_scl, Tongue_model, subjectMeta, featureMeta, model_path + "Tongue");

%%
% Plot PARAFAC models
%days = [-14 0 2 5 9 14 21];
days = 1:K;
timepoints = 1:K;
path_start = "./test_run/Figures/";

plotPARAFAC4(df_cnt_scl, Tongue_model, tongueVarExps{numFactors_tongue}, tongue_choice, subjectMeta, featureMeta, days, timepoints, 2, "PARAFAC tongue", path_start + "PARAFAC_tongue.jpg");

%%
% Show PCA performance
[coeff, score, latent, tsquared, explained] = pca(df_cnt_scl_flat, 'Centered',false, "NumComponents",numComponents);

subplot(4,2,1); bar(score(:,1)); % scores
subplot(4,2,2); bar(coeff(:,1)); % loadings
subplot(4,2,3); bar(score(:,2)); % scores
subplot(4,2,4); bar(coeff(:,2)); % loadings

subplot(4,2,5); bar(subjectLoadings(:,1));
subplot(4,2,6); bar(featureLoadings(:,1));
subplot(4,2,7); bar(subjectLoadings(:,2));
subplot(4,2,8); bar(featureLoadings(:,2));

%%
% Show parafac performance
[A,B,C] = fac2let(Tongue_model);
subplot(4,3,1); bar(B(:,1));
subplot(4,3,2); bar(A(:,1));
subplot(4,3,3); plot(1:K, C(:,1));
subplot(4,3,4); bar(B(:,2));
subplot(4,3,5); bar(A(:,2));
subplot(4,3,6); plot(1:K, C(:,2));

subplot(4,3,7); bar(featureLoadings(:,1));
subplot(4,3,8); bar(subjectLoadings(:,1));
subplot(4,3,9); plot(1:K, timeLoadings(:,1));
subplot(4,3,10); bar(featureLoadings(:,2));
subplot(4,3,11); bar(subjectLoadings(:,2));
subplot(4,3,12); plot(1:K, timeLoadings(:,2));