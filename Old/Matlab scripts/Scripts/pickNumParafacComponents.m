function numComponents=pickNumParafacComponents(X, maxFactors, numModelsPerComp, numRepsPerModel, maxCCdrop, minCCallowed)
Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 2;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)

core_consistency_mean = 1:maxFactors;
core_consistency_std = 1:maxFactors;

for i=1:maxFactors
    dummy = 1:numModelsPerComp;
    for j=1:numModelsPerComp
        [~, ~, ~, c] = myParafac(X, i, numRepsPerModel, Options);
        dummy(j) = c;
    end

    core_consistency_mean(i) = mean(dummy);
end

for k=2:maxFactors
    currentCoreConsistency = core_consistency_mean(k);
    prevCoreConsistency = core_consistency_mean(k-1);

    if(prevCoreConsistency-currentCoreConsistency>=maxCCdrop || currentCoreConsistency < minCCallowed)
        result = k;
        break;
    end
end

numComponents=result;