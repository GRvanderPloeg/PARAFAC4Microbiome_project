function [newAs, newBs, newCs]=checkComponentNegativity(realA, realB, realC, As, Bs, Cs)
I = size(As, 1);
J = size(Bs, 1);
K = size(Cs, 1);
numFactors = size(As, 2);
numReps = size(As, 3);

newAs = zeros(size(As));
newBs = zeros(size(Bs));
newCs = zeros(size(Cs));

for n=1:numFactors
    flatAs = reshape(As(:,n,:), I, numReps);
    flatBs = reshape(Bs(:,n,:), J, numReps);
    flatCs = reshape(Cs(:,n,:), K, numReps);

    meanA = realA(:,n);
    meanB = realB(:,n);
    meanC = realC(:,n);

    for r=1:numReps
        negEvidence = false(3,1);

        vectorA = flatAs(:,r);
        vectorB = flatBs(:,r);
        vectorC = flatCs(:,r);

        differenceA = sum((meanA - vectorA) .^ 2, "omitnan");
        negDifferenceA = sum((meanA + vectorA) .^ 2, "omitnan");
        differenceB = sum((meanB - vectorB) .^ 2, "omitnan");
        negDifferenceB = sum((meanB + vectorB) .^ 2, "omitnan");
        differenceC = sum((meanC - vectorC) .^ 2, "omitnan");
        negDifferenceC = sum((meanC + vectorC) .^ 2, "omitnan");

        if negDifferenceA < differenceA
            negEvidence(1) = true;
        end
        if negDifferenceB < differenceB
            negEvidence(2) = true;
        end
        if negDifferenceC < differenceC
            negEvidence(3) = true;
        end

        if sum(negEvidence) == 2
            if negEvidence(1) == 1
                newAs(:,n,r) = -1*vectorA;
            else
                newAs(:,n,r) = vectorA;
            end
            if negEvidence(2) == 1
                newBs(:,n,r) = -1*vectorB;
            else
                newBs(:,n,r) = vectorB;
            end
            if negEvidence(3) == 1
                newCs(:,n,r) = -1*vectorC;
            else
                newCs(:,n,r) = vectorC;
            end
        else
            newAs(:,n,r) = vectorA;
            newBs(:,n,r) = vectorB;
            newCs(:,n,r) = vectorC;
        end
    end
end

