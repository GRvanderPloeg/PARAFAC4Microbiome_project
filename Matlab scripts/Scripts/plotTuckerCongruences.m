function plotTuckerCongruences(allTuckers, overallTitle, path)

if ~exist('path','var')
    path="";
end

nModels = size(allTuckers{1},1);
nComparisons = size(allTuckers{1},2);

% Find the original number of components used by solving an equation
nFactors1 = (1 + sqrt(1+4*1*2*nComparisons)) / 2; % Quadratic equation formula
nFactors2 = (1 - sqrt(1+4*1*2*nComparisons)) / 2; % Quadratic equation formula

if nFactors1 > 0
    nComponents = nFactors1;
elseif nFactors2 > 0
    nComponents = nFactors2;
end

% Generate titles for all tucker comparisons made
comparisonTitles = strings(1,nComparisons);
comparisonIterator = 1;
for f1=1:(nComponents-1)
    for f2=(f1+1):nComponents
        comparisonTitles(comparisonIterator) = "Component " + f1 + " vs. " + f2;
        comparisonIterator = comparisonIterator + 1;
    end
end

% Make the plots
plotIterator = 1;

set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])

for j=1:nComparisons
    subplot(nComparisons,3,plotIterator); hold on;
    bar(1:nModels, allTuckers{1}(:,j));
    title("Tucker congruences of A (" + comparisonTitles(j) + ")");
    xlabel("Model number");
    ylabel("Tucker congruence coefficient");
    plotIterator = plotIterator + 1;
    
    subplot(nComparisons,3,plotIterator); hold on;
    bar(1:nModels, allTuckers{2}(:,j));
    title("Tucker congruences of B (" + comparisonTitles(j) + ")");
    xlabel("Model number");
    ylabel("Tucker congruence coefficient");
    plotIterator = plotIterator + 1;
    
    subplot(nComparisons,3,plotIterator); hold on;
    bar(1:nModels, allTuckers{3}(:,j));
    title("Tucker congruences of C (" + comparisonTitles(j) + ")");
    xlabel("Model number");
    ylabel("Tucker congruence coefficient");
    plotIterator = plotIterator + 1;
end

sgtitle(overallTitle);

if path ~= ""
    saveas(gcf, path)
    close()
end