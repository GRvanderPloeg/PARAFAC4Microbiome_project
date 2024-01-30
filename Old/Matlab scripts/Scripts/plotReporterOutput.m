function plotReporterOutput(allIterations, allErrors, allCorcondias, allVarExps, overallTitle, path)

if ~exist('path','var')
    path="";
end

nModels = size(allIterations,2);
set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])

subplot(2,2,1); hold on;
bar(1:nModels, allIterations);
ylim([0 max(allIterations)*1.2]);
title("Number of iterations needed per model");

subplot(2,2,2); hold on;
bar(1:nModels, allErrors);
ylim([0.9*min(allErrors) max(allErrors)*1.2]);
title("SSQ errors per model");

subplot(2,2,3); hold on;
bar(1:nModels, allCorcondias);
title("Corcondia per model");

subplot(2,2,4); hold on;
bar(1:nModels, allVarExps);
ylim([0 100]);
title("Variance explained per model");

sgtitle(overallTitle);

if path ~= ""
    saveas(gcf, path)
    close()
end