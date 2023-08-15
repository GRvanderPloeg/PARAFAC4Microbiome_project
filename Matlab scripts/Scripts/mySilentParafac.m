function [Factors,it,err,corcondia,stability,varExp,meanVarExp,stdVarExp]=mySilentParafac(X, numFactors, numReps, Options)
bestModel = [];
bestError = 1e50;
bestVarExp = 0;
allA = [];
allB = [];
allC = [];
allVarExp = 1:numReps;
allErrors = 1:numReps;
allowance = 0.0001;

for i=1:numReps
    %[Factors, it, err, c] = silent_parafac(X, numFactors, Options);
    evalc('[Factors, it, err, c] = silent_parafac(X, numFactors, Options)');
    [A, B, C] = fac2let(Factors);
    [A, B, C] = sortParafacComponents(X,A,B,C);

    allA(:,:,i) = A;
    allB(:,:,i) = B;
    allC(:,:,i) = C;
    allVarExp(i) = calcVarExplained(X, A, B, C);

    if err<bestError
        bestModel = Factors;
        bestIt = it;
        bestError = err;
        bestC = c;
        bestVarExp = calcVarExplained(X, A, B, C);
    end
end

[Abest,Bbest,Cbest] = fac2let(bestModel);
[Abest,Bbest,Cbest] = sortParafacComponents(X,Abest,Bbest,Cbest);
bestMatrix = Abest*krb(Cbest,Bbest)';

for k=1:numReps
    % Old stability metric: calculate difference between each vector in
    % each mode
    %Adiff = abs(Abest)/sum(abs(Abest)) - abs(allA(:,:,k))/sum(abs(allA(:,:,k)));
    %Bdiff = abs(Bbest) - abs(allB(:,:,k));
    %Cdiff = abs(Cbest) - abs(allC(:,:,k));
    %error = sum(Adiff) + sum(Bdiff) + sum(Cdiff);
    %error = sum(error);

    % Old stability metric: check the mean correlation on the diagonal of
    % the correlation matrix
    %M = allA(:,:,k)*krb(allC(:,:,k),allB(:,:,k))';
    %error = 1 - mean(diag(corr(bestMatrix, M)));
    %allErrors(k) = error;

    % Current stability metric: calculate modified RV for modeled matrices
    M = allA(:,:,k)*krb(allC(:,:,k),allB(:,:,k))';
    error = 1 - RV_modified(bestMatrix, M);
    allErrors(k) = error;
end

Factors = bestModel;
it = bestIt;
err = bestError;
corcondia = bestC;
stability = sum(allErrors <= allowance)/numReps*100;
varExp = bestVarExp;
meanVarExp = mean(allVarExp);
stdVarExp = std(allVarExp);