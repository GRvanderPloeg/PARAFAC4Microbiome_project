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

WA = (sampleInfo(:,7) == "WA" & sampleInfo(:,5) == "7");
WB = (sampleInfo(:,7) == "WB" & sampleInfo(:,5) == "8");
WC = (sampleInfo(:,7) == "WC" & sampleInfo(:,5) == "8");
SA = (sampleInfo(:,7) == "SA" & sampleInfo(:,5) == "7");
SB = (sampleInfo(:,7) == "SB" & sampleInfo(:,5) == "8");
SC = (sampleInfo(:,7) == "SC" & sampleInfo(:,5) == "4");

sampleSelection = WA | WB | WC | SA | SB | SC;
sampleInfo_filtered = sampleInfo(sampleSelection,:);

%%
% Filter based on sparsity per group

selection_WA = (sum(df(WA,:) == 0) / sum(WA));
selection_WB = (sum(df(WB,:) == 0) / sum(WB));
selection_WC = (sum(df(WC,:) == 0) / sum(WC));
selection_SA = (sum(df(SA,:) == 0) / sum(SA));
selection_SB = (sum(df(SB,:) == 0) / sum(SB));
selection_SC = (sum(df(SC,:) == 0) / sum(SC));

subplot(2,3,1); histogram(selection_WA);
subplot(2,3,2); histogram(selection_WB);
subplot(2,3,3); histogram(selection_WC);
subplot(2,3,4); histogram(selection_SA);
subplot(2,3,5); histogram(selection_SB);
subplot(2,3,6); histogram(selection_SC);

threshold = 0.75;
selection_WA = selection_WA <= threshold;
selection_WB = selection_WB <= threshold;
selection_WC = selection_WC <= threshold;
selection_SA = selection_SA <= threshold;
selection_SB = selection_SB <= threshold;
selection_SC = selection_SC <= threshold;

featureSelection = (selection_WA | selection_WB | selection_WC | selection_SA | selection_SB | selection_SC);
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
[df_cube, subjectMeta, conditionMeta] = rawDataToCube(df_clr, sampleInfo_filtered(:,7), str2double(sampleInfo_filtered(:,2)), keepIndividuals);

%%
% Center and scale
df_cnt = centerData(df_cube, 1);
df_cnt_scl = scaleData(df_cnt, 2);

%%
% Prepare metadata and remove outliers
subjectMeta(:,2) = ["Soil", "Soil", "Soil", "Water", "Water", "Water"];
subjectMeta(:,3) = ["Medium-A", "Medium-B", "Medium-C", "Medium-A", "Medium-B", "Medium-C"];

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
resort = [false false false];
legendIndex = [2 3 0];
numColsPerLegend = [2 3 0];
xlabels = ["Subject index" "Feature index" "Time index"];
titles = ["Subject mode", "Feature mode", "Time mode"];

[Models, Cons, VarExps, Boots, BootVarExps, Tuckers] = quickReport(df_cnt_scl_filtered, Options, const, maxComponents, numReps, maxIterations, metaData, resort, legendIndex, numColsPerLegend, xlabels, titles, "Fujita bootstrapped", path_start+"Fujita");

%%
% Dump
path_start = "./test_run/Dump/";
dump(Models, Cons, VarExps, Boots, BootVarExps, Tuckers, path_start, "Fujita");

%%
% Pick the number of components that seems appropriate based on the plots
% and pick the model that you like the best.
% Best practices: pick between best Corcondia or best variance explained.
% Alternative: choose a good tucker congruence model.

numFactors = 1;
choice = find(VarExps{numFactors}==max(VarExps{numFactors}));
Model = pickModel(Models, numFactors, choice);

%%
% Save models
model_path = "./test_run/PARAFAC models/";

annotatedModel = annotateModel(df_cnt_scl_filtered, Model, metaData);
savePARAFAC(df_cnt_scl_filtered, Model, annotatedModel, model_path + "SoH_microbiome");

%%
% Plot PARAFAC model
resort = [true true false];
legendIndex = [3 4 0];
numColsPerLegend = [2 3 0];
xlabels = ["Subject index" "Feature index" "Time index"];
titles = ["Subject mode", "Feature mode", "Time mode"];
path_start = "./test_run/Figures/";

plotPARAFAC4(annotatedModel, numFactors, VarExps, choice, resort, legendIndex, numColsPerLegend, xlabels, titles, "PARAFAC SoH microbiome", path_start + "PARAFAC_Fujita.jpg");
