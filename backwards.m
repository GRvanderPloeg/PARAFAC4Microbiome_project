% backwards
% Start from the very end, a finished PARAFAC model
addpath(".\Matlab scripts\Scripts\"); % own scripts
addpath(".\Matlab scripts\N-way toolbox\"); % from Rasmus Bro
%%
% Load loadings
A = readmatrix("./TIFN/Tongue_individual_mode.csv", Filetype="delimitedtext", OutputType="string");
B = readmatrix("./TIFN/Tongue_feature_mode.csv", Filetype="delimitedtext", OutputType="string");
C = readmatrix("./TIFN/Tongue_time_mode.csv", Filetype="delimitedtext", OutputType="string");
M = readmatrix("./TIFN/Tongue_model.csv", Filetype="delimitedtext");
Mcube = reshape(M, 41, 78, 7);
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

RESIDUALS = microb_tongue_cnt_scl - Mcube;

%%
% Initialize vectors
fakeA = str2double(A(:,1:2)); %normrnd(0, 11, [41 2]);
fakeB = str2double(B(:,1:2)); %normrnd(0.075, 0.075, [78 2]);
fakeC = str2double(C);

fakeM = fakeA * krb(fakeC, fakeB)';
fakeMcube = reshape(fakeM, [41 78 7]);
fakeMcube = fakeMcube + RESIDUALS;

%%
% Reverse scaling
fakeStds = microb_tongue_stds; %normrnd(1.9, 0.57, [78, 1]);
fakeMcube_revScl = fakeMcube;

for j=1:78
    for k=1:7
        fakeMcube_revScl(:,j,k) = fakeMcube_revScl(:,j,k) * fakeStds(j);
    end
end

%%
% Reverse centering
fakeMeans = microb_tongue_means; %normrnd(4, 1.72, [78, 7]);
fakeMcube_revScl_revCnt = fakeMcube_revScl;

for k=1:7
    for j=1:78
        fakeMcube_revScl_revCnt(:,j,k) = fakeMcube_revScl_revCnt(:,j,k) + fakeMeans(j,k);
    end
end

%%
% Reverse CLR
fakeGeoMeans = tongueGeoMeans; %normrnd(1.0328, 0.0082, 287);
dummy = permute(fakeMcube_revScl_revCnt, [1 3 2]);
fakeMfinal = reshape(dummy, 41*7, 78);

for i=1:287
    fakeMfinal(i,:) = exp(fakeMfinal(i,:)) * fakeGeoMeans(i);
end

fakeMfinal = round(fakeMfinal) - 1; % make into integers and remove pseudocount

%%
% Forward

% Replaced CLR as the procedure would not match the original one exactly
fakeMforward_clr = fakeMfinal+1;
for i=1:287
    fakeMforward_clr(i,:) = log(fakeMforward_clr(i,:) / tongueGeoMeans(i));
end

% Replace cube creation as the rows don't need to be reordered
fakeMforward_cube = reshape(fakeMforward_clr, 41, 78, 7);

% Rest is intact but not completely correct yet
[fakeMforward_cnt, fakeMforward_means] = centerData(fakeMforward_cube, 1);
[fakeMforward_cnt_scl, fakeMforward_stds] = scaleData(fakeMforward_cnt, 2);
pfac = parafac(fakeMforward_cnt_scl, 2);