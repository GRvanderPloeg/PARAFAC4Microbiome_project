function [allModels,allCons,allVarExps, bootstrappedModels, bootVarExps, allTuckers]=quickReport_SoH(X, maxComponents, numReps, individual_meta, timepoints, featureOffset, overallTitle, path_start)

allModels = cell(3,maxComponents);
allCons = cell(1,maxComponents);
allVarExps = cell(1,maxComponents);
bootstrappedModels = cell(3, maxComponents);
bootVarExps = cell(1,maxComponents);
allTuckers = cell(1, maxComponents);

for n=1:maxComponents
    [As, Bs, Cs, Its, Errs, Cons, varExps, Tuckers] = myReportingParafac(X, n, numReps);
    plotReporterOutput(Its, Errs, Cons, varExps, overallTitle, path_start+"_model_report_comp"+n+".jpg");
    if n >= 2
        plotTuckerCongruences(Tuckers, overallTitle, path_start+"_model_TCC_report_comp"+n+".jpg");
        allTuckers{n} = Tuckers;
    end
    allModels{1,n} = As;
    allModels{2,n} = Bs;
    allModels{3,n} = Cs;
    allCons{n} = Cons;
    allVarExps{n} = varExps;

    [As,Bs,Cs,varExps] = bootstrappedPARAFAC_SoH(X, individual_meta, n);
    bootstrappedModels{1,n} = As;
    bootstrappedModels{2,n} = Bs;
    bootstrappedModels{3,n} = Cs;
    bootVarExps{n} = varExps;
    plotBootstrappedPARAFAC_SoH(As, Bs, Cs, bootVarExps{n}, individual_meta, timepoints, featureOffset, overallTitle, path_start+"_bootstrapped_model_comp"+n+".jpg");
end