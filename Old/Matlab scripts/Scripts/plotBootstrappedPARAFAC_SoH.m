function plotBootstrappedPARAFAC_SoH(As, Bs, Cs, varExps, individual_mode_metadata, visitXlabel, featureOffset, title, path)

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

numFactors = size(As,2);

set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])

for i=0:(numFactors-1)

    Adata = reshape(As(:, i+1, :), size(As,1), size(As,3));
    individual_mode = [mean(Adata,2, 'omitnan') std(Adata,0,2, 'omitnan') individual_mode_metadata];
    individual_mode = sortrows(individual_mode, [4 3], "MissingPlacement","last");
    individual_mode = individual_mode(:,[1:2, size(individual_mode,2)]);
    individual_mode = [individual_mode (1:size(individual_mode,1))'];
    individual_mode = str2double(individual_mode);

    Bdata = reshape(Bs(:, i+1, :), size(Bs,1), size(Bs,3));
    feature_mode = [mean(Bdata,2) std(Bdata,0,2)];
    %feature_mode = str2double(feature_mode);
    Bdata
    feature_mode

    Cdata = reshape(Cs(:, i+1, :), size(Cs,1), size(Cs,3));
    time_mode = [mean(Cdata,2), std(Cdata,0,2)];

    % Feature vector
    subplot(numFactors, 3, 1+i*3), hold on;
    minY = min(feature_mode(:,1)-feature_mode(:,2));
    maxY = max(feature_mode(:,1)+feature_mode(:,2));
    temp = bar(feature_mode(:,1),"FaceColor","flat"); ylim([minY maxY]); ylabel("Comp. " + string(i+1)); xlabel("Feature index");

    % Color the bars
%     for j=1:max(feature_mode(:,3))
%         FeatureGroup = feature_mode(feature_mode(:,3)==j,:);
%         temp.CData(FeatureGroup(:,size(FeatureGroup,2)),:) = reshape(repelem(FeatureColors(j,:), size(FeatureGroup,1)), size(FeatureGroup,1), 3);
%     end

    % Add error bars
    er = errorbar(1:size(feature_mode,1), feature_mode(:,1), feature_mode(:,2), feature_mode(:,2));
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  

%     % Feature vector legend
%     h = zeros(max(feature_mode(:,3)),1);
%     for feature=1:max(feature_mode(:,3))
%         h(feature) = plot(NaN, NaN, "s", "MarkerEdgeColor", FeatureColors(feature,:), "MarkerFaceColor", FeatureColors(feature,:));
%     end
%     legend(h, feature_names, "Location", 'southoutside', "Orientation", "horizontal", "NumColumns", 3);
%     hold off;

    % Individual vector
    subplot(numFactors, 3, 2+i*3), hold on;
    minY = min(individual_mode(:,1) - individual_mode(:,2));
    maxY = max(individual_mode(:,1) + individual_mode(:,2));
    temp2 = bar(individual_mode(:,1),"FaceColor","flat"); ylim([minY maxY]); xlabel("Individual index");

    for k=0:max(individual_mode(:,3))
        RFgroup = individual_mode(individual_mode(:,3)==k,:);
        temp2.CData(RFgroup(:,size(RFgroup,2)),:) = reshape(repelem(RFcolors(k+1,:), size(RFgroup,1)), size(RFgroup,1), 3);
    end

    % Add error bars
    er = errorbar(1:size(individual_mode,1), individual_mode(:,1), individual_mode(:,2), individual_mode(:,2));
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  

    % ID vector legend
    h = zeros(3, 1);
    h(1) = plot(NaN,NaN,'s', "MarkerEdgeColor", RFcolors(1,:), "MarkerFaceColor", RFcolors(1,:));
    h(2) = plot(NaN,NaN,'s', "MarkerEdgeColor", RFcolors(2,:), "MarkerFaceColor", RFcolors(2,:));
    h(3) = plot(NaN,NaN,'s', "MarkerEdgeColor", RFcolors(3,:), "MarkerFaceColor", RFcolors(3,:));
    legend(h, {'Low response','Mid response','High response'}, "Location", 'southoutside', "Orientation", "horizontal");
    hold off;

    % Time vector
    subplot(numFactors, 3, 3+i*3), hold on;
    minY = min(time_mode(:,1) - time_mode(:,2));
    maxY = max(time_mode(:,1) + time_mode(:,2));
    plot(visitXlabel, time_mode(:,1)); ylim([minY maxY]); xlabel("Days"); xlim([min(visitXlabel) max(visitXlabel)]);
    yl = ylim;
    patch([0 14 14 0], [yl(1) yl(1) yl(2) yl(2)], 'red', "FaceAlpha", 0.3);
    errorbar(visitXlabel, time_mode(:,1), time_mode(:,2), "-", "Color", "black"); %xlabel("Days"); ylim([yl(1) yl(2)]);
    hold off;
    
end

sgtitle(title + ", variance explained = " + mean(varExps) + " +/- " + std(varExps) + " %");
saveas(gcf, path)
close()
end