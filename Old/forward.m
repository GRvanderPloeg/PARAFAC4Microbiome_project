% Forward
addpath(".\Matlab scripts\Scripts\"); % own scripts
addpath(".\Matlab scripts\N-way toolbox\"); % from Rasmus Bro

%%
% Load raw microbiome data
microbiome_raw = readmatrix("./TIFN/count-table.tsv", Filetype="delimitedtext", Delimiter="\t", OutputType="string");
microbiome_raw_controlgroup = microbiome_raw(microbiome_raw(:,4) == "control", :);

taxonomy = readmatrix("./TIFN/taxonomic-classification.tsv", Filetype="delimitedtext", Delimiter="\t", OutputType="string");
taxonomyHeader = taxonomy(1,:);
taxonomy = taxonomy(2:end,:);

subjectsControl = sortrows(unique(microbiome_raw(:,[2 4]), "rows"), 1);
subjectsControl = subjectsControl(subjectsControl(:,2) == "control",1);

%%
% Import red fluorescence data
rf_data = readmatrix("./TIFN/RFdata.csv", OutputType="string");
rf_data = rf_data(:, [1 6]);     % keep subject + RF group information

% Fix incorrect subject names
rf_data(rf_data(:,1) == "VSTPHZ", 1) = "VSTPH2";
rf_data(rf_data(:,1) == "D2VZH0", 1) = "DZVZH0";
rf_data(rf_data(:,1) == "DLODNN", 1) = "DLODDN";
rf_data(rf_data(:,1) == "O3VQFX", 1) = "O3VQFQ";
rf_data(rf_data(:,1) == "F80LGT", 1) = "F80LGF";
rf_data(rf_data(:,1) == "26QQR0", 1) = "26QQrO";

rf_data = unique(rf_data, "rows");
rf_data_control = rf_data(ismember(rf_data(:,1), subjectsControl), :);
rf_low = rf_data(rf_data(:,2) == "0", 1);
rf_mid = rf_data(rf_data(:,2) == "1", 1);
rf_high = rf_data(rf_data(:,2) == "2", 1);

%%
% Process data
microb_tongue_raw = microbiome_raw_controlgroup(microbiome_raw_controlgroup(:,5) == "tongue", :);
sparsityThreshold = 50;
[tongue_sparsity_low, tongue_sparsity_mid, tongue_sparsity_high] = calculateRFsparsity(microb_tongue_raw, rf_low, rf_mid, rf_high);
tongue_ASV_selection = (tongue_sparsity_low <= sparsityThreshold) | (tongue_sparsity_mid <= sparsityThreshold) | (tongue_sparsity_high <= sparsityThreshold);
nonbacterial = ~((taxonomy(:,5) == "Chloroplast") | (taxonomy(:,6) == "Mitochondria"));
microb_tongue_meta = microb_tongue_raw(:,1:5);
microb_tongue_numeric_strings = microb_tongue_raw(:,6:end);
microb_tongue_numeric = str2double(microb_tongue_numeric_strings);
[microb_tongue_clr, tongueGeoMeans] = transformCLR(microb_tongue_numeric);
microb_tongue_reduced = microb_tongue_clr(:, (tongue_ASV_selection & nonbacterial));
taxonomy_tongue_reduced = taxonomy((tongue_ASV_selection & nonbacterial), :);
numTimepoints_microb = 7;
microb_tongue = rawDataToCube_keepIndividuals(microb_tongue_reduced, microb_tongue_meta(:,2), microb_tongue_meta(:,3), numTimepoints_microb);
[microb_tongue_cnt, microb_tongue_means] = centerData(microb_tongue, 1);
[microb_tongue_cnt_scl, microb_tongue_stds] = scaleData(microb_tongue_cnt, 2);
microb_tongue_id_meta = rf_data_control(:,1);
microb_tongue_ASV_meta = taxonomy_tongue_reduced(:, [2:end 1]);

%%
% Determine array sizes
I = 41;
J = 78;
K = 9;
numComponents = 2;

%%
% Load fake data and process it
fakeMfinal = readmatrix("./temp.csv", FileType="delimitedtext");

% Replaced CLR as the procedure would not match the original one exactly
% Old way, using known geoMeans from the real data

% fakeMforward_clr = fakeMfinal+1;
% for i=1:287
%     fakeMforward_clr(i,:) = log(fakeMforward_clr(i,:) / tongueGeoMeans(i));
% end

% Proper way
[fakeMforward_clr, ~] = transformCLR(fakeMfinal);

% Replace cube creation as the rows don't need to be reordered
fakeMforward_cube = reshape(fakeMforward_clr, I, K, J);
fakeMforward_cube = permute(fakeMforward_cube, [1 3 2]);

% Rest is intact
[fakeMforward_cnt, fakeMforward_means] = centerData(fakeMforward_cube, 1);
[fakeMforward_cnt_scl, fakeMforward_stds] = scaleData(fakeMforward_cnt, 2);

%temp = microb_tongue_numeric(:, (tongue_ASV_selection & nonbacterial));

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

[tongueModels, tongueCons, tongueVarExps, tongueBoots, tongueBootVarExps, tongueTuckers] = quickReport(fakeMforward_cnt_scl, maxComponents, numReps, maxIterations, rf_data_control, microb_tongue_ASV_meta, days, 2, "Tongue bootstrapped", path_start+"tongue");

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
% Save the models
model_path = "./test_run/PARAFAC models/";

savePARAFAC(fakeMforward_cnt_scl, Tongue_model, microb_tongue_id_meta, microb_tongue_ASV_meta, model_path +  "Tongue");

%%
% Plot PARAFAC models
%days = [-14 0 2 5 9 14 21];
days = 1:K;
timepoints = 1:K;
path_start = "./test_run/Figures/";

plotPARAFAC4(fakeMforward_cnt_scl, Tongue_model, tongueVarExps{numFactors_tongue}, tongue_choice, rf_data_control, microb_tongue_ASV_meta, days, timepoints, 2, "PARAFAC tongue", path_start + "PARAFAC_tongue.jpg");
