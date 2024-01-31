% PARAFAC functionality
addpath("..\N-way-shell\Scripts\"); % own scripts
addpath("..\N-way-shell\N-way toolbox/"); % from Rasmus Bro

%%
% Load data
df = readmatrix("./Data/16S_data.csv");
taxonomy = readmatrix("./Data/taxonomy.csv", OutputType="string");
taxonomy = taxonomy(2:end,:);
sampleInfo = readmatrix("./Data/sampleInfo.csv", OutputType="string");

%%
% Select only specific replicates
WC = (sampleInfo(:,7) == "WC");
sampleSelection = WC;
sampleInfo_filtered = sampleInfo(WC,:);

%%
% Filter based on sparsity per group
selection_WC = (sum(df(WC,:) == 0) / sum(WC));
threshold = 0.99;
featureSelection = selection_WC < threshold;
taxonomy_filtered = taxonomy(featureSelection,:);

%%
% Filter df
df_filtered = df(sampleSelection,featureSelection);

%%
% CLR
df_clr = transformCLR(df_filtered);

%%
% Make into cube
keepIndividuals = true;
[df_cube, subjectMeta, conditionMeta] = rawDataToCube(df_clr, str2double(sampleInfo_filtered(:,5)), str2double(sampleInfo_filtered(:,2)), keepIndividuals);

%%
% Center and scale
df_cnt = centerData(df_cube, 1);
df_cnt_scl = scaleData(df_cnt, 2);

%%
% Prepare metadata and remove outliers
subjectMeta = [subjectMeta subjectMeta];
featureMeta = taxonomy_filtered;
conditionMeta = conditionMeta;

% Remove outliers here
subjectMeta_filtered = subjectMeta;
featureMeta_filtered = featureMeta;
conditionMeta_filtered = conditionMeta;
df_cnt_scl_filtered = df_cnt_scl;

%%
% Initialize options
Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 2;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)
const = [0 0 0];

%%
% Run PARAFACs
path_start = "./test_run/Figures/";
maxComponents=5;
numReps=25;
maxIterations = 20;
metaData = {subjectMeta_filtered, featureMeta_filtered, conditionMeta_filtered};
resort = [false true false];
legendIndex = [0 3 0];
numColsPerLegend = [0 3 0];
xlabels = ["Subject index" "Feature index" "Time index"];
titles = ["Subject mode", "Feature mode", "Time mode"];
balancedJackKnifing = false;
subjectGroupCol = 0;

[Models, Cons, VarExps, Boots, BootVarExps, Tuckers] = quickReport(df_cnt_scl_filtered, Options, const, maxComponents, numReps, maxIterations, balancedJackKnifing, subjectGroupCol, metaData, resort, legendIndex, numColsPerLegend, xlabels, titles, "Fujita bootstrapped", path_start+"Fujita");

%%
% Dump
path_start = "./test_run/Dump/";
dump(Models, Cons, VarExps, Boots, BootVarExps, Tuckers, path_start, "Fujita");

%%
% Pick the number of components that seems appropriate based on the plots
% and pick the model that you like the best.
% Best practices: pick between best Corcondia or best variance explained.
% Alternative: choose a good tucker congruence model.

numFactors = 3;
choice = find(VarExps{numFactors}==max(VarExps{numFactors}));
Model = pickModel(Models, numFactors, choice);

%%
% Save models
model_path = "./test_run/PARAFAC models/";

annotatedModel = annotateModel(df_cnt_scl_filtered, Model, metaData);
savePARAFAC(df_cnt_scl_filtered, Model, annotatedModel, model_path + "SoH_microbiome");

%%
% Plot PARAFAC model
resort = [false true false];
legendIndex = [0 6 0];
numColsPerLegend = [0 3 0];
xlabels = ["Subject index" "Feature index" "Time index"];
titles = ["Subject mode", "Feature mode", "Time mode"];
path_start = "./test_run/Figures/";

plotPARAFAC4(annotatedModel, numFactors, VarExps, choice, resort, legendIndex, numColsPerLegend, xlabels, titles, "PARAFAC Fujita", path_start + "PARAFAC_Fujita.jpg");
