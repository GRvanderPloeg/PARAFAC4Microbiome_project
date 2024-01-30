function permutedPARAFAC(X, numComponents, numReps, path)
numModelsToGenerate_real = 10; % generate x models for a good distr
numModelsToGenerate_perm = 1; % generate x models for a quick varExp

Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 2;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)

numIDs = size(X,1);
numFeatures = size(X,2);
numTimepoints = size(X,3);

% Normal execution of PARAFAC
[~,~,~,~,~,~,meanVarExp_real,stdVarExp_real] = mySilentParafac(X, numComponents, numModelsToGenerate_real, Options);

% To generate variances explained for the real model without permuting, we
% will use the variance explained distribution of the real model.
gm = gmdistribution(meanVarExp_real, stdVarExp_real);
varExps_real = random(gm, numReps);

disp("real modeling complete!")

% Permutation of individuals per feature, time mode kept intact
% (aka shifting tubes around per vertical slice).
varExps_id_perm = 1:numReps;

for i=1:numReps
    X_id_permuted = X;
    for j=1:numFeatures
        ordering = randperm(numIDs);
        for k=1:numTimepoints
            X_id_permuted(:,j,k) = X_id_permuted(ordering,j,k);
        end
    end
    
    [~,~,~,~,~,varExp_id_perm,~,~] = mySilentParafac(X_id_permuted, numComponents, numModelsToGenerate_perm, Options);
    varExps_id_perm(i) = varExp_id_perm;
end

disp("id permutation modeling complete!")

% Permutation of individuals per feature and per timepoint
% aka breaking everything
varExps_id_time_perm = 1:numReps;

for i=1:numReps
    X_id_time_permuted = X;
    for j=1:numFeatures
        for k=1:numTimepoints
            ordering = randperm(numIDs);
            X_id_time_permuted(:,j,k) = X_id_time_permuted(ordering,j,k);
        end
    end
    
    [~,~,~,~,~,varExp_id_time_perm,~,~] = mySilentParafac(X_id_time_permuted, numModelsToGenerate_perm, 10, Options);
    varExps_id_time_perm(i) = varExp_id_time_perm;
end

disp("id and time permutation modeling complete!")

% Save the variance explained distributions
writematrix(varExps_real, path+"_varExps_real.csv");
writematrix(varExps_id_perm', path+"_varExps_id_perm.csv");
writematrix(varExps_id_time_perm', path+"_varExps_id_time_perm.csv");