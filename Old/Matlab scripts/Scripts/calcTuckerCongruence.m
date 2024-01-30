function [tuckersA, tuckersB, tuckersC] = calcTuckerCongruence(A, B, C)
numFactors = size(A,2);

if numFactors == 1
    %disp('You cannot calculate similarity in a 1-component model.');
    tuckersA = 0;
    tuckersB = 0;
    tuckersC = 0;
    return;
end

numCombinations = numFactors * (numFactors-1) / 2;
numCombo = 1;
tuckersA = zeros(1,numCombinations);
tuckersB = zeros(1,numCombinations);
tuckersC = zeros(1,numCombinations);

for f1=1:(numFactors-1)
    for f2=(f1+1):numFactors
        tuckersA(numCombo) = sum(A(:,f1) .* A(:,f2)) / sqrt(sumsqr(A(:,f1)) * sumsqr(A(:,f2)));
        tuckersB(numCombo) = sum(B(:,f1) .* B(:,f2)) / sqrt(sumsqr(B(:,f1)) * sumsqr(B(:,f2)));
        tuckersC(numCombo) = sum(C(:,f1) .* C(:,f2)) / sqrt(sumsqr(C(:,f1)) * sumsqr(C(:,f2)));
        numCombo = numCombo + 1;
    end
end