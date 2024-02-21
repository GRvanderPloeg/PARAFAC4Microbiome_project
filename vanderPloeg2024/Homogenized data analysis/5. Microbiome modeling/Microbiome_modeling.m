% PARAFAC functionality
addpath("..\..\..\N-way-shell\Scripts\"); % own scripts
addpath("..\..\..\N-way-shell\N-way toolbox\"); % from Rasmus Bro

%%
% Load input data
numTimepoints_microb = 7;

microb_tongue_cnt_scl = readmatrix("../2. Microbiome pre-processing/20230605_run/Microb_tongue.csv");
microb_tongue_cnt_scl = reshape(microb_tongue_cnt_scl, size(microb_tongue_cnt_scl, 1), size(microb_tongue_cnt_scl,2)/numTimepoints_microb, numTimepoints_microb);

microb_lowling_cnt_scl = readmatrix("../2. Microbiome pre-processing/20230605_run/Microb_lowling.csv");
microb_lowling_cnt_scl = reshape(microb_lowling_cnt_scl, size(microb_lowling_cnt_scl, 1), size(microb_lowling_cnt_scl,2)/numTimepoints_microb, numTimepoints_microb);

microb_lowinter_cnt_scl = readmatrix("../2. Microbiome pre-processing/20230605_run/Microb_lowinter.csv");
microb_lowinter_cnt_scl = reshape(microb_lowinter_cnt_scl, size(microb_lowinter_cnt_scl, 1), size(microb_lowinter_cnt_scl,2)/numTimepoints_microb, numTimepoints_microb);

microb_upling_cnt_scl = readmatrix("../2. Microbiome pre-processing/20230605_run/Microb_upling.csv");
microb_upling_cnt_scl = reshape(microb_upling_cnt_scl, size(microb_upling_cnt_scl, 1), size(microb_upling_cnt_scl,2)/numTimepoints_microb, numTimepoints_microb);

microb_upinter_cnt_scl = readmatrix("../2. Microbiome pre-processing/20230605_run/Microb_upinter.csv");
microb_upinter_cnt_scl = reshape(microb_upinter_cnt_scl, size(microb_upinter_cnt_scl, 1), size(microb_upinter_cnt_scl,2)/numTimepoints_microb, numTimepoints_microb);

microb_saliva_cnt_scl = readmatrix("../2. Microbiome pre-processing/20230605_run/Microb_saliva.csv");
microb_saliva_cnt_scl = reshape(microb_saliva_cnt_scl, size(microb_saliva_cnt_scl, 1), size(microb_saliva_cnt_scl,2)/numTimepoints_microb, numTimepoints_microb);

%%
% Load metadata
microb_tongue_id_meta = readmatrix("../2. Microbiome pre-processing/20230605_run/Tongue_id_meta.csv", OutputType="string");
microb_lowling_id_meta = readmatrix("../2. Microbiome pre-processing/20230605_run/Lowling_id_meta.csv", OutputType="string");
microb_lowinter_id_meta = readmatrix("../2. Microbiome pre-processing/20230605_run/Lowinter_id_meta.csv", OutputType="string");
microb_upling_id_meta = readmatrix("../2. Microbiome pre-processing/20230605_run/Upling_id_meta.csv", OutputType="string");
microb_upinter_id_meta = readmatrix("../2. Microbiome pre-processing/20230605_run/Upinter_id_meta.csv", OutputType="string");
microb_saliva_id_meta = readmatrix("../2. Microbiome pre-processing/20230605_run/Saliva_id_meta.csv", OutputType="string");

microb_tongue_ASV_meta = readmatrix("../2. Microbiome pre-processing/20230605_run/Tongue_feature_meta.csv", OutputType="string");
microb_lowling_ASV_meta = readmatrix("../2. Microbiome pre-processing/20230605_run/Lowling_feature_meta.csv", OutputType="string");
microb_lowinter_ASV_meta = readmatrix("../2. Microbiome pre-processing/20230605_run/Lowinter_feature_meta.csv", OutputType="string");
microb_upling_ASV_meta = readmatrix("../2. Microbiome pre-processing/20230605_run/Upling_feature_meta.csv", OutputType="string");
microb_upinter_ASV_meta = readmatrix("../2. Microbiome pre-processing/20230605_run/Upinter_feature_meta.csv", OutputType="string");
microb_saliva_ASV_meta = readmatrix("../2. Microbiome pre-processing/20230605_run/Saliva_feature_meta.csv", OutputType="string");

%%
% Import red fluorescence data
rf_data = readmatrix("../0. Raw data input/RFdata.csv", OutputType="string");
rf_data = rf_data(:, [1 6]);     % keep subject + RF group information

% Fix incorrect subject names
rf_data(rf_data(:,1) == "VSTPHZ", 1) = "VSTPH2";
rf_data(rf_data(:,1) == "D2VZH0", 1) = "DZVZH0";
rf_data(rf_data(:,1) == "DLODNN", 1) = "DLODDN";
rf_data(rf_data(:,1) == "O3VQFX", 1) = "O3VQFQ";
rf_data(rf_data(:,1) == "F80LGT", 1) = "F80LGF";
rf_data(rf_data(:,1) == "26QQR0", 1) = "26QQrO";

rf_data = unique(rf_data, "rows");
rf_data = rf_data(ismember(rf_data(:,1), microb_tongue_id_meta),:);

%%
% Initialize options
Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 2;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)
const = [0 0 0];

%%
% Prepare metadata
rf_data_mod = rf_data;
rf_data_mod(rf_data(:,2) == "0", 2) = "Low";
rf_data_mod(rf_data(:,2) == "1", 2) = "Mid";
rf_data_mod(rf_data(:,2) == "2", 2) = "High";

metaData_tongue = {};
metaData_tongue{1} = rf_data_mod; % subjects
metaData_tongue{2} = microb_tongue_ASV_meta; % features
metaData_tongue{3} = ["Baseline" "Baseline", "Intervention", "Intervention", "Intervention", "Intervention", "Resolution"]';      % timepoints

metaData_lowling = {};
metaData_lowling{1} = rf_data_mod; % subjects
metaData_lowling{2} = microb_lowling_ASV_meta; % features
metaData_lowling{3} = ["Baseline" "Baseline", "Intervention", "Intervention", "Intervention", "Intervention", "Resolution"]';      % timepoints

metaData_lowinter = {};
metaData_lowinter{1} = rf_data_mod; % subjects
metaData_lowinter{2} = microb_lowinter_ASV_meta; % features
metaData_lowinter{3} = ["Baseline" "Baseline", "Intervention", "Intervention", "Intervention", "Intervention", "Resolution"]';      % timepoints

metaData_upling = {};
metaData_upling{1} = rf_data_mod; % subjects
metaData_upling{2} = microb_upling_ASV_meta; % features
metaData_upling{3} = ["Baseline" "Baseline", "Intervention", "Intervention", "Intervention", "Intervention", "Resolution"]';      % timepoints

metaData_upinter = {};
metaData_upinter{1} = rf_data_mod; % subjects
metaData_upinter{2} = microb_upinter_ASV_meta; % features
metaData_upinter{3} = ["Baseline" "Baseline", "Intervention", "Intervention", "Intervention", "Intervention", "Resolution"]';      % timepoints

metaData_saliva = {};
metaData_saliva{1} = rf_data_mod; % subjects
metaData_saliva{2} = microb_saliva_ASV_meta; % features
metaData_saliva{3} = ["Baseline" "Baseline", "Intervention", "Intervention", "Intervention", "Intervention", "Resolution"]';      % timepoints

%%
% Do a reporter Parafac + plot it to check what the best model is.
% path_start = "./20230904_run_newCodeBase/Figures/";
% maxComponents=3;
% days = [-14 0 2 5 9 14 21];
% numReps=25;
% maxIterations=20;
% 
% [tongueModels, tongueCons, tongueVarExps, tongueBoots, tongueBootVarExps, tongueTuckers] = quickReport(microb_tongue_cnt_scl, maxComponents, numReps, maxIterations, rf_data, microb_tongue_ASV_meta, days, 2, "Tongue bootstrapped", path_start+"tongue");
% [lowlingModels, lowlingCons, lowlingVarExps, lowlingBoots, lowlingBootVarExps, lowlingTuckers] = quickReport(microb_lowling_cnt_scl, maxComponents, numReps, maxIterations, rf_data, microb_lowling_ASV_meta, days, 2, "Low ling bootstrapped", path_start+"lowling");
% [lowinterModels, lowinterCons, lowinterVarExps, lowinterBoots, lowinterBootVarExps, lowinterTuckers] = quickReport(microb_lowinter_cnt_scl, maxComponents, numReps, maxIterations, rf_data, microb_lowinter_ASV_meta, days, 2, "Low inter bootstrapped", path_start+"lowinter");
% [uplingModels, uplingCons, uplingVarExps, uplingBoots, uplingBootVarExps, uplingTuckers] = quickReport(microb_upling_cnt_scl, maxComponents, numReps, maxIterations, rf_data, microb_upling_ASV_meta, days, 2, "Up ling bootstrapped", path_start+"upling");
% [upinterModels, upinterCons, upinterVarExps, upinterBoots, upinterBootVarExps, upinterTuckers] = quickReport(microb_upinter_cnt_scl, maxComponents, numReps, maxIterations, rf_data, microb_upinter_ASV_meta, days, 2, "Up inter bootstrapped", path_start+"upinter");
% [salivaModels, salivaCons, salivaVarExps, salivaBoots, salivaBootVarExps, salivaTuckers] = quickReport(microb_saliva_cnt_scl, maxComponents, numReps, maxIterations, rf_data, microb_saliva_ASV_meta, days, 2, "Saliva bootstrapped", path_start+"saliva");

path_start = "./test_run/Figures/";
maxComponents=4;
days = [-14 0 2 5 9 14 21];
numReps=50;
maxIterations=20;
resort = [true true false];
legendIndex = [2 2 1]; % which column of the annotatedModel mode has the metadata to colour by.
numColsPerLegend = [3 2 3];
xlabels = ["Subject index", "Feature index", "Time index"];
titles = ["Subject mode", "Feature mode", "Time mode"];
balancedJackKnifing = true;
subjectGroupCol = 2;

[tongueModels, tongueCons, tongueVarExps, tongueBoots, tongueBootVarExps, tongueTuckers] = quickReport(microb_tongue_cnt_scl, Options, const, maxComponents, numReps, maxIterations, balancedJackKnifing, subjectGroupCol, metaData_tongue, resort, legendIndex, numColsPerLegend, xlabels, titles, "Tongue bootstrapped", path_start+"tongue");
[lowlingModels, lowlingCons, lowlingVarExps, lowlingBoots, lowlingBootVarExps, lowlingTuckers] = quickReport(microb_lowling_cnt_scl, Options, const, maxComponents, numReps, maxIterations, balancedJackKnifing, subjectGroupCol, metaData_lowling, resort, legendIndex, numColsPerLegend, xlabels, titles, "Low ling bootstrapped", path_start+"lowling");
[lowinterModels, lowinterCons, lowinterVarExps, lowinterBoots, lowinterBootVarExps, lowinterTuckers] = quickReport(microb_lowinter_cnt_scl, Options, const, maxComponents, numReps, maxIterations, balancedJackKnifing, subjectGroupCol, metaData_lowinter, resort, legendIndex, numColsPerLegend, xlabels, titles, "Low inter bootstrapped", path_start+"lowinter");
[uplingModels, uplingCons, uplingVarExps, uplingBoots, uplingBootVarExps, uplingTuckers] = quickReport(microb_upling_cnt_scl, Options, const, maxComponents, numReps, maxIterations, balancedJackKnifing, subjectGroupCol, metaData_upling, resort, legendIndex, numColsPerLegend, xlabels, titles, "Up ling bootstrapped", path_start+"upling");
[upinterModels, upinterCons, upinterVarExps, upinterBoots, upinterBootVarExps, upinterTuckers] = quickReport(microb_upinter_cnt_scl, Options, const, maxComponents, numReps, maxIterations, balancedJackKnifing, subjectGroupCol, metaData_upinter, resort, legendIndex, numColsPerLegend, xlabels, titles, "Up inter bootstrapped", path_start+"upinter");
[salivaModels, salivaCons, salivaVarExps, salivaBoots, salivaBootVarExps, salivaTuckers] = quickReport(microb_saliva_cnt_scl, Options, const, maxComponents, numReps, maxIterations, balancedJackKnifing, subjectGroupCol, metaData_saliva, resort, legendIndex, numColsPerLegend, xlabels, titles, "Saliva bootstrapped", path_start+"saliva");

%%
% Dump data so far for later inspection
path_start = "./test_run/Dump/";

dump(tongueModels, tongueCons, tongueVarExps, tongueBoots, tongueBootVarExps, tongueTuckers, path_start, "tongue");
dump(lowlingModels, lowlingCons, lowlingVarExps, lowlingBoots, lowlingBootVarExps, lowlingTuckers, path_start, "lowling");
dump(lowinterModels, lowinterCons, lowinterVarExps, lowinterBoots, lowinterBootVarExps, lowinterTuckers, path_start, "lowinter");
dump(uplingModels, uplingCons, uplingVarExps, uplingBoots, uplingBootVarExps, uplingTuckers, path_start, "upling");
dump(upinterModels, upinterCons, upinterVarExps, upinterBoots, upinterBootVarExps, upinterTuckers, path_start, "upinter");
dump(salivaModels, salivaCons, salivaVarExps, salivaBoots, salivaBootVarExps, salivaTuckers, path_start, "saliva");

%%
% Pick the number of components that seems appropriate based on the plots
% and pick the model that you like the best.
% Best practices: pick between best Corcondia or best variance explained.
% Alternative: choose a good tucker congruence model.

numFactors_tongue = 2;
numFactors_lowling = 2;
numFactors_lowinter = 1;
numFactors_upling = 2;
numFactors_upinter = 2;
numFactors_saliva = 2;

tongue_choice = find(tongueVarExps{numFactors_tongue}==max(tongueVarExps{numFactors_tongue}));
%tongue_choice = find(tongueCons{numFactors_tongue}==max(tongueCons{numFactors_tongue}));

lowling_choice = find(lowlingVarExps{numFactors_lowling}==max(lowlingVarExps{numFactors_lowling}));
%lowling_choice = find(lowlingCons{numFactors_lowling}==max(lowlingCons{numFactors_lowling}));

lowinter_choice = find(lowinterVarExps{numFactors_lowinter}==max(lowinterVarExps{numFactors_lowinter}));
%lowinter_choice = find(lowinterCons{numFactors_lowinter}==max(lowinterCons{numFactors_lowinter}));

upling_choice = find(uplingVarExps{numFactors_upling}==max(uplingVarExps{numFactors_upling}));
%upling_choice = find(uplingCons{numFactors_upling}==max(uplingCons{numFactors_upling}));

upinter_choice = find(upinterVarExps{numFactors_upinter}==max(upinterVarExps{numFactors_upinter}));
%upinter_choice = find(upinterCons{numFactors_upinter}==max(upinterCons{numFactors_upinter}));

saliva_choice = find(salivaVarExps{numFactors_saliva}==max(salivaVarExps{numFactors_saliva}));

Tongue_model = pickModel(tongueModels, numFactors_tongue, tongue_choice);
L_ling_model = pickModel(lowlingModels, numFactors_lowling, lowling_choice);
L_inter_model = pickModel(lowinterModels, numFactors_lowinter, lowinter_choice);
U_ling_model = pickModel(uplingModels, numFactors_upling, upling_choice);
U_inter_model = pickModel(upinterModels, numFactors_upinter, upinter_choice);
Saliva_model = pickModel(salivaModels, numFactors_saliva, saliva_choice);

%%
% Save the models
% model_path = "./20230904_run_newCodeBase/PARAFAC models/";
% 
% savePARAFAC(microb_tongue_cnt_scl, Tongue_model, microb_tongue_id_meta, microb_tongue_ASV_meta, model_path +  "Tongue");
% savePARAFAC(microb_lowinter_cnt_scl, L_inter_model, microb_lowinter_id_meta, microb_lowinter_ASV_meta, model_path +  "Low_inter");
% savePARAFAC(microb_lowling_cnt_scl, L_ling_model, microb_lowling_id_meta, microb_lowling_ASV_meta, model_path + "Low_ling");
% savePARAFAC(microb_upinter_cnt_scl, U_inter_model, microb_upinter_id_meta, microb_upinter_ASV_meta, model_path + "Up_inter");
% savePARAFAC(microb_upling_cnt_scl, U_ling_model, microb_upling_id_meta, microb_upling_ASV_meta, model_path + "Up_ling");
% savePARAFAC(microb_saliva_cnt_scl, Saliva_model, microb_saliva_id_meta, microb_saliva_ASV_meta, model_path + "Saliva");

model_path = "./test_run/PARAFAC models/";

annotatedModel_tongue = annotateModel(microb_tongue_cnt_scl, Tongue_model, metaData_tongue);
annotatedModel_lowling = annotateModel(microb_lowling_cnt_scl, L_ling_model, metaData_lowling);
annotatedModel_lowinter = annotateModel(microb_lowinter_cnt_scl, L_inter_model, metaData_lowinter);
annotatedModel_upling = annotateModel(microb_upling_cnt_scl, U_ling_model, metaData_upling);
annotatedModel_upinter = annotateModel(microb_upinter_cnt_scl, U_inter_model, metaData_upinter);
annotatedModel_saliva = annotateModel(microb_saliva_cnt_scl, Saliva_model, metaData_saliva);

savePARAFAC(microb_tongue_cnt_scl, Tongue_model, annotatedModel_tongue, model_path + "Tongue");
savePARAFAC(microb_lowling_cnt_scl, L_ling_model, annotatedModel_lowling, model_path + "Low_ling");
savePARAFAC(microb_lowinter_cnt_scl, L_inter_model, annotatedModel_lowinter, model_path + "Low_inter");
savePARAFAC(microb_upling_cnt_scl, U_ling_model, annotatedModel_upling, model_path + "Up_ling");
savePARAFAC(microb_upinter_cnt_scl, U_inter_model, annotatedModel_upinter, model_path + "Up_inter");
savePARAFAC(microb_saliva_cnt_scl, Saliva_model, annotatedModel_saliva, model_path + "Saliva");

%%
% Plot PARAFAC models
% days = [-14 0 2 5 9 14 21];
% timepoints = 1:7;
% path_start = "./20230904_run_newCodeBase/Figures/";
% 
% plotPARAFAC4(microb_tongue_cnt_scl, Tongue_model, tongueVarExps{numFactors_tongue}, tongue_choice, rf_data, microb_tongue_ASV_meta, days, timepoints, 2, "PARAFAC tongue", path_start + "PARAFAC_tongue.jpg");
% plotPARAFAC4(microb_lowinter_cnt_scl, L_inter_model, lowinterVarExps{numFactors_lowling}, lowinter_choice, rf_data, microb_lowinter_ASV_meta, days, timepoints, 2, "PARAFAC lower jaw, interproximal", path_start + "PARAFAC_low_inter.jpg");
% plotPARAFAC4(microb_lowling_cnt_scl, L_ling_model, lowlingVarExps{numFactors_lowinter}, lowling_choice, rf_data, microb_lowling_ASV_meta, days, timepoints, 2, "PARAFAC lower jaw, lingual", path_start + "PARAFAC_low_ling.jpg");
% plotPARAFAC4(microb_upinter_cnt_scl, U_inter_model, upinterVarExps{numFactors_upinter}, upinter_choice, rf_data, microb_upinter_ASV_meta, days, timepoints, 2, "PARAFAC upper jaw, interproximal", path_start + "PARAFAC_up_inter.jpg");
% plotPARAFAC4(microb_upling_cnt_scl, U_ling_model, uplingVarExps{numFactors_upling}, upling_choice, rf_data, microb_upling_ASV_meta, days, timepoints, 2, "PARAFAC upper jaw, lingual", path_start + "PARAFAC_up_ling.jpg");
% plotPARAFAC4(microb_saliva_cnt_scl, Saliva_model, salivaVarExps{numFactors_saliva}, saliva_choice, rf_data, microb_saliva_ASV_meta, days, timepoints, 2, "PARAFAC saliva", path_start + "PARAFAC_saliva.jpg");

legendIndex2 = [4 4 3]; % which column of the annotatedModel mode has the metadata to colour by.
legendIndex1 = [3 3 2];
numColsPerLegend = [3 2 3];
resort = [true true false];
xlabels = ["Subject index", "Feature index", "Time index"];
titles = ["Subject mode", "Feature mode", "Time mode"];
path_start = "./test_run/Figures/";

plotPARAFAC4(annotatedModel_tongue, numFactors_tongue, tongueVarExps, tongue_choice, resort, legendIndex2, numColsPerLegend, xlabels, titles, "PARAFAC tongue", path_start + "PARAFAC_tongue.jpg");
plotPARAFAC4(annotatedModel_lowling, numFactors_lowling, lowlingVarExps, lowling_choice, resort, legendIndex2, numColsPerLegend, xlabels, titles, "PARAFAC lower jaw, lingual", path_start + "PARAFAC_low_ling.jpg");
plotPARAFAC4(annotatedModel_lowinter, numFactors_lowinter, lowinterVarExps, lowinter_choice, resort, legendIndex1, numColsPerLegend, xlabels, titles, "PARAFAC lower jaw, interproximal", path_start + "PARAFAC_low_inter.jpg");
plotPARAFAC4(annotatedModel_upling, numFactors_upling, uplingVarExps, upling_choice, resort, legendIndex2, numColsPerLegend, xlabels, titles, "PARAFAC upper jaw, lingual", path_start + "PARAFAC_up_ling.jpg");
plotPARAFAC4(annotatedModel_upinter, numFactors_upinter, upinterVarExps, upinter_choice, resort, legendIndex2, numColsPerLegend, xlabels, titles, "PARAFAC upper jaw, interproximal", path_start + "PARAFAC_up_inter.jpg");
plotPARAFAC4(annotatedModel_saliva, numFactors_saliva, salivaVarExps, saliva_choice, resort, legendIndex2, numColsPerLegend, xlabels, titles, "PARAFAC saliva", path_start + "PARAFAC_saliva.jpg");