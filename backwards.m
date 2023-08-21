% backwards
% Start from the very end, a finished PARAFAC model
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
% Prepare loadings
A = readmatrix("./TIFN/Tongue_individual_mode.csv", Filetype="delimitedtext", OutputType="string");
B = readmatrix("./TIFN/Tongue_feature_mode.csv", Filetype="delimitedtext", OutputType="string");
C = readmatrix("./TIFN/Tongue_time_mode.csv", Filetype="delimitedtext", OutputType="string");

% Convoluted way of generating similar-ish subject loadings to tongue data.
realAnum = str2double(A(:,[1 2]));
Anum = realAnum;
%Anum = normrnd(0, 11, [41 2]);

Ameans = [1 0.5 -3; -1 -0.5 3]; % rows is component, cols is RFgroup
Astds = [10 10 10; 10 10 10];   % rows is component, cols is RFgroup

mask_low = rf_data_control(:,2) == "0";
mask_mid = rf_data_control(:,2) == "1";
mask_high = rf_data_control(:,2) == "2";

% Anum(mask_low, 1) = normrnd(mean(realAnum(mask_low, 1)), std(realAnum(mask_low, 1)), size(realAnum(mask_low,1)));
% Anum(mask_low, 2) = normrnd(mean(realAnum(mask_low, 2)), std(realAnum(mask_low, 2)), size(realAnum(mask_low,2)));
% Anum(mask_mid, 1) = normrnd(mean(realAnum(mask_mid, 1)), std(realAnum(mask_mid, 1)), size(realAnum(mask_mid,1)));
% Anum(mask_mid, 2) = normrnd(mean(realAnum(mask_mid, 2)), std(realAnum(mask_mid, 2)), size(realAnum(mask_mid,2)));
% Anum(mask_high, 1) = normrnd(mean(realAnum(mask_high, 1)), std(realAnum(mask_high, 1)), size(realAnum(mask_high,1)));
% Anum(mask_high, 2) = normrnd(mean(realAnum(mask_high, 2)), std(realAnum(mask_high, 2)), size(realAnum(mask_high,2)));

Anum(mask_low, 1) = normrnd(Ameans(1,1), Astds(1,1), size(realAnum(mask_low,1)));
Anum(mask_low, 2) = normrnd(Ameans(2,1), Astds(2,1), size(realAnum(mask_low,2)));
Anum(mask_mid, 1) = normrnd(Ameans(1,2), Astds(1,2), size(realAnum(mask_mid,1)));
Anum(mask_mid, 2) = normrnd(Ameans(2,2), Astds(2,2), size(realAnum(mask_mid,2)));
Anum(mask_high, 1) = normrnd(Ameans(1,3), Astds(1,3), size(realAnum(mask_high,1)));
Anum(mask_high, 2) = normrnd(Ameans(2,3), Astds(2,3), size(realAnum(mask_high,2)));

% Convoluted way to generating similar-ish feature loadings to tongue data.
realBnum = str2double(B(:,[1 2]));
Bnum = realBnum;
%Bnum = normrnd(0.075, 0.075, [78 2]);

%phyla = unique(B(:,4));
mask_actinobacteriota = B(:,4) == "Actinobacteriota";
mask_bacteroidota = B(:,4) == "Bacteroidota";
mask_campilobacterota = B(:,4) == "Campilobacterota";
mask_firmicutes = B(:,4) == "Firmicutes";
mask_fusobacteriota = B(:,4) == "Fusobacteriota";
mask_patescibacteria = B(:,4) == "Patescibacteria";
mask_proteobacteria = B(:,4) == "Proteobacteria";

% Bnum(mask_actinobacteriota, 1) = normrnd(mean(realBnum(mask_actinobacteriota, 1)), std(realBnum(mask_actinobacteriota, 1)), size(realBnum(mask_actinobacteriota, 1)));
% Bnum(mask_actinobacteriota, 2) = normrnd(mean(realBnum(mask_actinobacteriota, 2)), std(realBnum(mask_actinobacteriota, 2)), size(realBnum(mask_actinobacteriota, 2)));
% Bnum(mask_bacteroidota, 1) = normrnd(mean(realBnum(mask_bacteroidota, 1)), std(realBnum(mask_bacteroidota, 1)), size(realBnum(mask_bacteroidota, 1)));
% Bnum(mask_bacteroidota, 2) = normrnd(mean(realBnum(mask_bacteroidota, 2)), std(realBnum(mask_bacteroidota, 2)), size(realBnum(mask_bacteroidota, 2)));
% Bnum(mask_campilobacterota, 1) = normrnd(mean(realBnum(mask_campilobacterota, 1)), std(realBnum(mask_campilobacterota, 1)), size(realBnum(mask_campilobacterota, 1)));
% Bnum(mask_campilobacterota, 2) = normrnd(mean(realBnum(mask_campilobacterota, 2)), std(realBnum(mask_campilobacterota, 2)), size(realBnum(mask_campilobacterota, 2)));
% Bnum(mask_firmicutes, 1) = normrnd(mean(realBnum(mask_firmicutes, 1)), std(realBnum(mask_firmicutes, 1)), size(realBnum(mask_firmicutes, 1)));
% Bnum(mask_firmicutes, 2) = normrnd(mean(realBnum(mask_firmicutes, 2)), std(realBnum(mask_firmicutes, 2)), size(realBnum(mask_firmicutes, 2)));
% Bnum(mask_fusobacteriota, 1) = normrnd(mean(realBnum(mask_fusobacteriota, 1)), std(realBnum(mask_fusobacteriota, 1)), size(realBnum(mask_fusobacteriota, 1)));
% Bnum(mask_fusobacteriota, 2) = normrnd(mean(realBnum(mask_fusobacteriota, 2)), std(realBnum(mask_fusobacteriota, 2)), size(realBnum(mask_fusobacteriota, 2)));
% Bnum(mask_patescibacteria, 1) = normrnd(mean(realBnum(mask_patescibacteria, 1)), std(realBnum(mask_patescibacteria, 1)), size(realBnum(mask_patescibacteria, 1)));
% Bnum(mask_patescibacteria, 2) = normrnd(mean(realBnum(mask_patescibacteria, 2)), std(realBnum(mask_patescibacteria, 2)), size(realBnum(mask_patescibacteria, 2)));
% Bnum(mask_proteobacteria, 1) = normrnd(mean(realBnum(mask_proteobacteria, 1)), std(realBnum(mask_proteobacteria, 1)), size(realBnum(mask_proteobacteria, 1)));
% Bnum(mask_proteobacteria, 2) = normrnd(mean(realBnum(mask_proteobacteria, 2)), std(realBnum(mask_proteobacteria, 2)), size(realBnum(mask_proteobacteria, 2)));

Bmeans = [0.05 0.1 -0.05 0.05 0.2 0.15 -0.1; 0.05 0.15 0.1 0.075 0.15 -0.1 0.1];
Bstds = [0.1 0.1 0 0.1 0.1 0 0.05; 0.1 0.05 0 0.1 0.05 0 0.05];

Bnum(mask_actinobacteriota, 1) = normrnd(Bmeans(1,1), Bstds(1,1), size(realBnum(mask_actinobacteriota, 1)));
Bnum(mask_actinobacteriota, 2) = normrnd(Bmeans(2,1), Bstds(2,1), size(realBnum(mask_actinobacteriota, 2)));
Bnum(mask_bacteroidota, 1) = normrnd(Bmeans(1,2), Bstds(1,2), size(realBnum(mask_bacteroidota, 1)));
Bnum(mask_bacteroidota, 2) = normrnd(Bmeans(2,2), Bstds(2,2), size(realBnum(mask_bacteroidota, 2)));
Bnum(mask_campilobacterota, 1) = normrnd(Bmeans(1,3), Bstds(1,3), size(realBnum(mask_campilobacterota, 1)));
Bnum(mask_campilobacterota, 2) = normrnd(Bmeans(2,3), Bstds(2,3), size(realBnum(mask_campilobacterota, 2)));
Bnum(mask_firmicutes, 1) = normrnd(Bmeans(1,4), Bstds(1,4), size(realBnum(mask_firmicutes, 1)));
Bnum(mask_firmicutes, 2) = normrnd(Bmeans(2,4), Bstds(2,4), size(realBnum(mask_firmicutes, 2)));
Bnum(mask_fusobacteriota, 1) = normrnd(Bmeans(1,5), Bstds(1,5), size(realBnum(mask_fusobacteriota, 1)));
Bnum(mask_fusobacteriota, 2) = normrnd(Bmeans(2,5), Bstds(2,5), size(realBnum(mask_fusobacteriota, 2)));
Bnum(mask_patescibacteria, 1) = normrnd(Bmeans(1,6), Bstds(1,6), size(realBnum(mask_patescibacteria, 1)));
Bnum(mask_patescibacteria, 2) = normrnd(Bmeans(2,6), Bstds(2,6), size(realBnum(mask_patescibacteria, 2)));
Bnum(mask_proteobacteria, 1) = normrnd(Bmeans(1,7), Bstds(1,7), size(realBnum(mask_proteobacteria, 1)));
Bnum(mask_proteobacteria, 2) = normrnd(Bmeans(2,7), Bstds(2,7), size(realBnum(mask_proteobacteria, 2)));

% time
%Cnum = str2double(C(:,[1 2]));
Cnum = [0.1 0.9; 0.11 0.8; 0.2 0.65; 0.55 0.55; 0.3 0.5; 0.25 0.44; 0.2 0.3; 0.11 0.25; 0.1 0.65];

%%
% Determine array sizes
I = size(Anum, 1);
J = size(Bnum, 1);
K = size(Cnum, 1);
numComponents = size(Anum, 2);

%%
% Create processed matrix

% Original modelled data matrix
%M = readmatrix("./TIFN/Tongue_component_1.csv", Filetype="delimitedtext");
%Mcube = reshape(M, 41, 78, 7);

% Calculate residuals
%RESIDUALS = microb_tongue_cnt_scl - Mcube;
%RESIDUALS = normrnd(0, 0.8265, [41 78 7]);

% Or generate from vectors
fakeM = Anum * krb(Cnum,Bnum)';
fakeMcube = reshape(fakeM, I, J, K);

% Add residual structure
%fakeMcube = fakeMcube + RESIDUALS;

%%
% Reverse scaling
%fakeStds = microb_tongue_stds; %normrnd(1.9, 0.57, [78, 1]);
fakeStds = normrnd(1.9, 0.57, [J 1]);

fakeMcube_revScl = fakeMcube;

for j=1:J
    for k=1:K
        fakeMcube_revScl(:,j,k) = fakeMcube_revScl(:,j,k) * fakeStds(j);
    end
end

%%
% Reverse centering
%fakeMeans = microb_tongue_means; %normrnd(4, 1.72, [78, 7]);
fakeMeans = normrnd(4, 1.72, [J K]);

fakeMcube_revScl_revCnt = fakeMcube_revScl;

for k=1:K
    for j=1:J
        fakeMcube_revScl_revCnt(:,j,k) = fakeMcube_revScl_revCnt(:,j,k) + fakeMeans(j,k);
    end
end

%%
% Reverse CLR
%fakeGeoMeans = tongueGeoMeans;
fakeGeoMeans = normrnd(1.0328, 0.0082, [I*K 1]);

dummy = permute(fakeMcube_revScl_revCnt, [1 3 2]);
fakeMfinal = reshape(dummy, I*K, J);

for i=1:(I*K)
    fakeMfinal(i,:) = exp(fakeMfinal(i,:)) * fakeGeoMeans(i);
end

fakeMfinal = round(fakeMfinal) - 1; % make into integers and remove pseudocount

%%
% Replace negative counts with zero
min(fakeMfinal, [], "all")
mask = fakeMfinal < 0;
fakeMfinal(mask) = 0;

%%
% Save
writematrix(fakeMfinal, "./temp.csv");