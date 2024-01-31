% PARAFAC functionality
addpath("..\N-way-shell\Scripts\"); % own scripts
addpath("..\N-way-shell\N-way toolbox/"); % from Rasmus Bro

%%
% Load data
df = readmatrix("./Data/count-table_clean.csv");
taxonomy = readmatrix("./Data/speciesMetadata_clean.csv", OutputType="string");
sampleInfo = readmatrix("./Data/sampleMetadata_clean.csv", OutputType="string");

% Filter samples to just the infants and the right timepoints
sampleMask1 = (sampleInfo(:,6) ~= "Mother"); % conserves late infant samples
sampleMask2 = ismember(sampleInfo(:,3), ["4", "7", "21", "Infancy"]);
sampleMask2 = ismember(sampleInfo(:,3), ["4", "7", "21"]);
sampleMask = sampleMask1 & sampleMask2;
sampleInfo = sampleInfo(sampleMask,:);
%sampleInfo(sampleInfo(:,3) == "Infancy", 3) = "24";

% Filter taxa to taxa with at least 1 non-zero value
featureMask = (sum(df) > 0);
df = df(sampleMask, featureMask);
taxonomy = taxonomy(featureMask,:);

%%
% Filter based on sparsity per group
threshold = 0.90;
vaginallyBorn = sampleInfo(:,4) == "Vaginal";
caesareanBorn = sampleInfo(:,4) == "Caesarean";

vaginalDF = df(vaginallyBorn,:);
caesareanDF = df(caesareanBorn,:);

vaginalSparsity = sum(vaginalDF == 0, 1) / size(vaginalDF, 1);
caesareanSparsity = sum(caesareanDF == 0, 1) / size(caesareanDF, 1);

%subplot(1,2,1); histogram(vaginalSparsity);
%subplot(1,2,2); histogram(caesareanSparsity);

vaginalSelection = vaginalSparsity <= threshold;
caesareanSelection = caesareanSparsity <= threshold;
featureSelection = vaginalSelection | caesareanSelection;

taxonomy_filtered = taxonomy(featureSelection,:);
df_filtered = df(:,featureSelection);

%%
% CLR
df_clr = transformCLR(df_filtered);

%%
% Make into cube
keepIndividuals = false;
[df_cube, subjectMeta, conditionMeta] = rawDataToCube(df_clr, sampleInfo(:,2), str2double(sampleInfo(:,3)), keepIndividuals);

%%
% Center and scale
df_cnt = centerData(df_cube, 1);
df_cnt_scl = scaleData(df_cnt, 2);

%%
% Prepare metadata and remove outliers
subjectMeta_filtered = unique(sampleInfo(:,[2 4]), "rows");
subjectMeta_filtered = subjectMeta_filtered(ismember(subjectMeta_filtered(:,1), subjectMeta),:);
featureMeta_filtered = taxonomy_filtered;
timeMeta_filtered = conditionMeta;

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
maxComponents=3;
numReps=25;
maxIterations = 20;
metaData = {subjectMeta_filtered, featureMeta_filtered, timeMeta_filtered};
resort = [true false false];
legendIndex = [2 0 0];
numColsPerLegend = [2 0 0];
xlabels = ["Subject index" "Feature index" "Time index"];
titles = ["Subject mode", "Feature mode", "Time mode"];
balancedJackKnifing = true;
subjectGroupCol = 2;

[Models, Cons, VarExps, Boots, BootVarExps, Tuckers] = quickReport(df_cnt_scl_filtered, Options, const, maxComponents, numReps, maxIterations, balancedJackKnifing, subjectGroupCol, metaData, resort, legendIndex, numColsPerLegend, xlabels, titles, "Shao bootstrapped", path_start+"Shao");

%%
% Dump
path_start = "./test_run/Dump/";
dump(Models, Cons, VarExps, Boots, BootVarExps, Tuckers, path_start, "Shao");

%%
% Pick the number of components that seems appropriate based on the plots
% and pick the model that you like the best.
% Best practices: pick between best Corcondia or best variance explained.
% Alternative: choose a good tucker congruence model.

numFactors = 2;
choice = find(Cons{numFactors}==max(Cons{numFactors}));
Model = pickModel(Models, numFactors, choice);

%%
% Save models
model_path = "./test_run/PARAFAC models/";

annotatedModel = annotateModel(df_cnt_scl_filtered, Model, metaData);
savePARAFAC(df_cnt_scl_filtered, Model, annotatedModel, model_path + "Shao");

%%
% Plot PARAFAC model
resort = [true false false];
legendIndex = [4 0 0];
numColsPerLegend = [2 0 0];
xlabels = ["Subject index" "Feature index" "Time index"];
titles = ["Subject mode", "Feature mode", "Time mode"];
path_start = "./test_run/Figures/";

plotPARAFAC4(annotatedModel, numFactors, VarExps, choice, resort, legendIndex, numColsPerLegend, xlabels, titles, "PARAFAC Shao", path_start + "PARAFAC_Shao.jpg");
