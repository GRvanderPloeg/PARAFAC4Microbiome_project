function [sortedA,sortedB,sortedC]=sortParafacComponents(X, A, B, C)
numFactors = size(A,2);
varExplained = 1:numFactors;

for i=1:numFactors
    varExplained(i) = calcVarExplained(X, A(:,i), B(:,i), C(:,i));
end

[~, Index] = sort(varExplained, 'descend');

sortedA = A(:,Index);
sortedB = B(:,Index);
sortedC = C(:,Index);