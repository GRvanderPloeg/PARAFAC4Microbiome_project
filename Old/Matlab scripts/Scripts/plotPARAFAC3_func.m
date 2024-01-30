function plotPARAFAC3_func(X, model, stability, varExps, choiceNumber, individual_mode_metadata, feature_mode_metadata, visitXlabel, usedVisits, title, path)

% FeatureOffset = 1 for metabolomics (if you want super_pathway colors)
% = 2 for microbiome (if you want phylum colors)

FeatureColors =    [122 62 118; % rgb(122, 62, 118)
                    61 105 123; % rgb(61, 105, 123)
                    49 155 49; % rgb(49, 155, 49)
                    183 107 1; % rgb(183, 107, 1)
                    184 3 0; % rgb(184, 3, 0)
                    180 106 175; % rgb(180, 106, 175)
                    91 151 174; % rgb(91, 151, 174)
                    84 201 84; % rgb(84, 201, 84)
                    254 153 11; % rgb(254, 153, 11)
                    255 34 31; % rgb(255, 34, 31)
                    204 153 201; % rgb(204, 153, 201)
                    158 193 207; % rgb(158, 193, 207)
                    158 224 158; % rgb(158, 224, 158)
                    254 177 68; % rgb(254, 177, 68)
                    255 102 99; % rgb(255, 102, 99)
                    255 0 0;
                    0 255 0;
                    0 0 255;
                    0 0 0;
                    0 255 255;
                    255 0 255;
                    128 128 128;
                    128 0 0;
                    128 128 0;
                    0 128 0;
                    128 0 128;
                    0 128 128;
                    0 0 128;
                    191 191 191;
                    0 0 0]/255;

RFcolors = [0.7529 0.7529 0.7529;
            0 0 1;
            1 0 0];

[A, B, C] = fac2let(model);
[A, B, C] = sortParafacComponents(X, A, B, C);
numFactors = size(A,2);

individual_mode = [A individual_mode_metadata];
individual_mode = sortrows(individual_mode, size(A,2)+2, "MissingPlacement","last");
individual_mode = individual_mode(:,[1:numFactors, size(individual_mode,2)]);
individual_mode = [individual_mode (1:size(individual_mode,1))'];
individual_mode = str2double(individual_mode);

feature_mode = [B feature_mode_metadata];
feature_mode = sortrows(feature_mode, (size(B,2)+1):size(feature_mode,2), "MissingPlacement","last");
feature_mode = [feature_mode(:,1:size(B,2)) feature_mode(:,size(B,2)+1)];
feature_mode = [feature_mode (1:size(feature_mode,1))'];
feature_names = unique(feature_mode(:,numFactors+1));
feature_mode(:,size(feature_mode,2)-1) = grp2idx(feature_mode(:,size(feature_mode,2)-1));
feature_mode = str2double(feature_mode);

set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])

for i=0:(numFactors-1)
    varExplained = calcVarExplained(X, A(:,i+1), B(:,i+1), C(:,i+1));

    % Feature vector
    subplot(numFactors, 3, 1+i*3), hold on;
    temp = bar(feature_mode(:,i+1),"FaceColor","flat"); ylabel("Comp. " + string(i+1) + " (" + string(varExplained) + "%)"); xlabel("Feature index");

    % Individual vector
    subplot(numFactors, 3, 2+i*3), hold on;
    temp2 = bar(individual_mode(:,i+1),"FaceColor","flat"); xlabel("Individual index");

    for k=0:max(individual_mode(:,size(individual_mode,2)-1))
        RFgroup = individual_mode(individual_mode(:,size(individual_mode,2)-1)==k,:);
        temp2.CData(RFgroup(:,size(RFgroup,2)),:) = reshape(repelem(RFcolors(k+1,:), size(RFgroup,1)), size(RFgroup,1), 3);
    end

    % ID vector legend
    h = zeros(3, 1);
    h(1) = plot(NaN,NaN,'s', "MarkerEdgeColor", RFcolors(1,:), "MarkerFaceColor", RFcolors(1,:));
    h(2) = plot(NaN,NaN,'s', "MarkerEdgeColor", RFcolors(2,:), "MarkerFaceColor", RFcolors(2,:));
    h(3) = plot(NaN,NaN,'s', "MarkerEdgeColor", RFcolors(3,:), "MarkerFaceColor", RFcolors(3,:));
    legend(h, {'Low response','Mid response','High response'}, "Location", 'southoutside', "Orientation", "horizontal");

    % Time vector
    subplot(numFactors, 3, 3+i*3), hold on;
    plot(visitXlabel, C(usedVisits,i+1)); xlabel("Days"); xlim([min(visitXlabel) max(visitXlabel)]);
    yl = ylim;
    patch([0 14 14 0], [yl(1) yl(1) yl(2) yl(2)], 'red', "FaceAlpha", 0.3);
    plot(visitXlabel, C(usedVisits,i+1), "-", "Color", "black"); xlabel("Days"); ylim([yl(1) yl(2)]);
    scatter(visitXlabel, C(usedVisits,i+1), "filled", "MarkerEdgeColor", "black", "MarkerFaceColor","black");

end

varExp = varExps(choiceNumber);
meanVarExp = mean(varExps);
stdVarExp = std(varExps);
sgtitle(title + ", stability=" + stability + "%, variance explained =  " + varExp + "% (" + meanVarExp + " +/- " + stdVarExp + "%)");
saveas(gcf, path)
close()
end