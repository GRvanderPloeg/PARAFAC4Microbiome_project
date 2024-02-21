% PARAFAC functionality
home = ".";
cd(home)
addpath("..\Matlab scripts\Scripts\"); % own scripts
addpath("..\Matlab scripts\N-way toolbox\"); % from Rasmus Bro

%%
numComponents = 2;
numTimepoints = 7;
inputPath = "../5. Microbiome modeling/20230618_run/PARAFAC models/Saliva";
path = "./20230618_run/Saliva";

input_data = str2double(readmatrix(inputPath + "_input.csv", Filetype="delimitedtext", Delimiter=",", OutputType="string"));
[I, J] = size(input_data);
J = J / numTimepoints;
K = numTimepoints;
input_data = reshape(input_data, I, J, K);

individual_mode = readmatrix(inputPath + "_individual_mode.csv", Filetype="delimitedtext", Delimiter=",", OutputType="string");
feature_mode = readmatrix(inputPath + "_feature_mode.csv", Filetype="delimitedtext", Delimiter=",", OutputType="string");
time_mode = readmatrix(inputPath + "_time_mode.csv", Filetype="delimitedtext", Delimiter=",", OutputType="string");

A = str2double(individual_mode(:,1:numComponents));
B = str2double(feature_mode(:,1:numComponents));
C = str2double(time_mode(:,1:numComponents));
pfac_model = {A, B, C};

congruence_result_individuals = conload(input_data, pfac_model, 1);
congruence_result_features = conload(input_data, pfac_model, 2);
congruence_result_time = conload(input_data, pfac_model, 3);

writematrix(congruence_result_individuals, path + "_individual_congruence_loadings.csv");
writematrix(congruence_result_features, path + "_feature_congruence_loadings.csv");
writematrix(congruence_result_time, path + "_time_congruence_loadings.csv");
