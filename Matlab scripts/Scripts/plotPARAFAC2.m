function plotPARAFAC2(X, numFactors, numReps, individual_mode_metadata, feature_mode_metadata, visitXlabel, usedVisits, title, path)

Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 2;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)

FeatureColors = [0 0.4470 0.7410;
                0.8500 0.3250 0.0980;
                0.9290 0.6940 0.1250;
                0.4940 0.1840 0.5560;
                0.4660 0.6740 0.1880;
                0.3010 0.7450 0.9330;
                0.6350 0.0780 0.1840;
                1 1 1;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0];

RFcolors = [0.7529 0.7529 0.7529;
            0 0 1;
            1 0 0];

[Factors, ~, ~, ~, stability, varExp, meanVarExp, stdVarExp] = myParafac(X, numFactors, numReps, Options);
[A, B, C] = fac2let(Factors);
[A, B, C] = sortParafacComponents(X, A, B, C);

individual_mode = [A individual_mode_metadata];
individual_mode = sortrows(individual_mode, (size(A,2)+1):size(individual_mode,2), "MissingPlacement","last");
individual_mode = individual_mode(:,1:(size(individual_mode,2)-1));
individual_mode = [individual_mode (1:size(individual_mode,1))'];
individual_mode = str2double(individual_mode);

feature_mode = [B feature_mode_metadata];
feature_mode = sortrows(feature_mode, (size(B,2)+1):size(feature_mode,2), "MissingPlacement","last");
feature_mode = [feature_mode(:,1:size(B,2)) feature_mode(:,size(B,2)+1)];
feature_mode = [feature_mode (1:size(feature_mode,1))'];
feature_mode(:,size(feature_mode,2)-1) = grp2idx(feature_mode(:,size(feature_mode,2)-1));
feature_mode = str2double(feature_mode);

set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])

for i=0:(numFactors-1)
    varExplained = calcVarExplained(X, A(:,i+1), B(:,i+1), C(:,i+1));

    % Feature vector
    subplot(numFactors, 3, 1+i*3), hold on;
    temp = bar(feature_mode(:,i+1),"FaceColor","flat"); ylabel("Comp. " + string(i+1) + " (" + string(varExplained) + "%)");

    for j=1:max(feature_mode(:,size(feature_mode,2)-1))
        FeatureGroup = feature_mode(feature_mode(:,(size(feature_mode,2)-1))==j,:);
        temp.CData(FeatureGroup(:,size(FeatureGroup,2)),:) = reshape(repelem(FeatureColors(j,:), size(FeatureGroup,1)), size(FeatureGroup,1), 3);
    end

    % Individual vector
    subplot(numFactors, 3, 2+i*3), hold on;
    temp2 = bar(individual_mode(:,i+1),"FaceColor","flat");

    for k=1:max(individual_mode(:,size(individual_mode,2)-1))
        RFgroup = individual_mode(individual_mode(:,size(individual_mode,2)-1)==k,:);
        temp2.CData(RFgroup(:,size(RFgroup,2)),:) = reshape(repelem(RFcolors(k,:), size(RFgroup,1)), size(RFgroup,1), 3);
    end

    % Time vector
    subplot(numFactors, 3, 3+i*3), plot(visitXlabel, C(usedVisits,i+1));
end

sgtitle(title + ", stability=" + stability + "%, variance explained =  " + varExp + "% (" + meanVarExp + " +/- " + stdVarExp + "%)");
saveas(gcf, path)
close()
end