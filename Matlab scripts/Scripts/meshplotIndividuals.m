function meshplotIndividuals(X, numFactors, numReps, individual_metadata, path)

Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 2;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)

model = myParafac(X,numFactors,numReps,Options);
M = nmodel(model);
E = X - M;

LowRF = datasample(find(str2double(individual_metadata(:,1)) == 1), 3, "Replace", false);
MidRF = datasample(find(str2double(individual_metadata(:,1)) == 2), 3, "Replace", false);
HighRF = datasample(find(str2double(individual_metadata(:,1)) == 3), 3, "Replace", false);
individuals = [LowRF' MidRF' HighRF'];

plotDim = ceil(sqrt(size(individuals,2)));
%set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8.12, 14.44], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])
set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])

for i=1:size(individuals,2)
    id = individuals(i);
    subplot(plotDim,plotDim,i),mesh(reshape(X(id,:,:),[size(X,2) size(X,3)]))
    xlabel("Timepoints");
    ylabel("Features");
    zlabel("Values");

    if i <= 3
        title("Low");
    elseif i <= 6
        title("Mid");
    else
        title("High");
    end
end
saveas(gcf,path+"_raw.jpg");
close()

set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])

for i=1:size(individuals,2)
    id = individuals(i);
    subplot(plotDim,plotDim,i),mesh(reshape(E(id,:,:),[size(X,2) size(X,3)]))
    xlabel("Timepoints");
    ylabel("Features");
    zlabel("Values");

    if i <= 3
        title("Low");
    elseif i <= 6
        title("Mid");
    else
        title("High");
    end
end
saveas(gcf,path+"_residuals.jpg");
close()

set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])

for i=1:size(individuals,2)
    id = individuals(i);
    subplot(plotDim,plotDim,i),mesh(reshape(M(id,:,:),[size(X,2) size(X,3)]))
    xlabel("Timepoints");
    ylabel("Features");
    zlabel("Values");

    if i <= 3
        title("Low");
    elseif i <= 6
        title("Mid");
    else
        title("High");
    end
end
saveas(gcf,path+"_model.jpg");
close()