function [allA,allB,allC,allIterations,allErrors,allCorcondias,allVarExp,allTuckers]=myReportingParafac(X, numFactors, numReps)
Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 2;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)
Options(4) = 2;     % no scaling (default 1)
Options(5) = NaN;	% no output (default 10)
Options(6) = 2500;   % max number of iterations (default 2500)
const = [0 0 0];    % PARAFAC model constraints, currently none

allA = [];
allB = [];
allC = [];
allIterations = 1:numReps;
allErrors = 1:numReps;
allCorcondias = 1:numReps;
allVarExp = 1:numReps;
allTuckers = cell(1,3);

textprogressbar('Creating PARAFAC models ');
for i=1:numReps
    textprogressbar(i/numReps*100);
    %[Factors, it, err, c] = silent_parafac(X, numFactors, Options);
    evalc('[Factors, it, err, c] = silent_parafac(X, numFactors, Options, const)');
    [A, B, C] = fac2let(Factors);

    allA(:,:,i) = A;
    allB(:,:,i) = B;
    allC(:,:,i) = C;
    allIterations(i) = it;
    allErrors(i) = err;
    allCorcondias(i) = c;
    allVarExp(i) = calcVarExplained(X, A, B, C);
    [tuckersA, tuckersB, tuckersC] = calcTuckerCongruence(A, B, C);
    allTuckers{1}(i,:) = tuckersA;
    allTuckers{2}(i,:) = tuckersB;
    allTuckers{3}(i,:) = tuckersC;
end

textprogressbar(' done');