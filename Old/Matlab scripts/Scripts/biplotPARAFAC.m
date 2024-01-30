function biplotPARAFAC(X, model, individual_mode_metadata, feature_mode_metadata, featureOffset, plot_title, path)
% Only supports 2 and 3-component PARAFAC models!

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
feature_mode = [feature_mode(:,1:size(B,2)) feature_mode(:,size(B,2)+featureOffset)];
feature_mode = [feature_mode (1:size(feature_mode,1))'];

if numFactors==2
    feature_names = unique(feature_mode(:,3));
elseif numFactors==3
    feature_names = unique(feature_mode(:,4));
end

feature_mode(:,size(feature_mode,2)-1) = grp2idx(feature_mode(:,size(feature_mode,2)-1));
feature_mode = str2double(feature_mode);

set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])

if numFactors==2

    % Feature vector
    subplot(1,3,1),hold on;
    xlabel("Feature vector 1");
    ylabel("Feature vector 2");
    Vector1 = feature_mode(:,1);
    Vector2 = feature_mode(:,2);
    angle = acos((Vector1'*Vector2)/(norm(Vector1)*norm(Vector2)));
    title("Feature mode, angle = " + (angle/(2*pi))*360 + " degrees");
    for i=1:max(feature_mode(:,3))
        plot(feature_mode(feature_mode(:,3)==i,1), feature_mode(feature_mode(:,3)==i,2), ".", "MarkerSize", 20, "Color", FeatureColors(i,:));
    end

    % feature legend
    h = zeros(max(feature_mode(:,numFactors+1)),1);
    for feature=1:max(feature_mode(:,numFactors+1))
        h(feature) = plot(NaN, NaN, "s", "MarkerEdgeColor", FeatureColors(feature,:), "MarkerFaceColor", FeatureColors(feature,:));
    end
    legend(h, feature_names, "Location", 'southoutside', "Orientation", "horizontal", "NumColumns", 3);

    % Individual vectors
    subplot(1,3,2),hold on;
    xlabel("Individual vector 1");
    ylabel("Individual vector 2");
    Vector1 = individual_mode(:,1);
    Vector2 = individual_mode(:,2);
    angle = acos((Vector1'*Vector2)/(norm(Vector1)*norm(Vector2)));
    title("Individual mode, angle = " + (angle/(2*pi))*360 + " degrees");
    for i=0:2
        plot(individual_mode(individual_mode(:,3) == i,1),individual_mode(individual_mode(:,3) == i,2),".", "MarkerSize", 20, "Color",RFcolors(i+1,:));
    end

    % Add id legend
    h = zeros(3, 1);
    h(1) = plot(NaN,NaN,'s', "MarkerEdgeColor", RFcolors(1,:), "MarkerFaceColor", RFcolors(1,:));
    h(2) = plot(NaN,NaN,'s', "MarkerEdgeColor", RFcolors(2,:), "MarkerFaceColor", RFcolors(2,:));
    h(3) = plot(NaN,NaN,'s', "MarkerEdgeColor", RFcolors(3,:), "MarkerFaceColor", RFcolors(3,:));
    legend(h, {'Low response','Mid response','High response'}, "Location", 'southoutside', "Orientation", "horizontal");
    
    % Time vectors
    labels = ["-14", "0", "2", "5", "9", "14", "21"];
    labels2 = ["0", "2", "5", "9", "14"];

    subplot(1,3,3),hold on;
    xlabel("Time vector 1");
    ylabel("Time vector 2");
    Vector1 = C(:,1);
    Vector2 = C(:,2);
    angle = acos((Vector1'*Vector2)/(norm(Vector1)*norm(Vector2)));
    title("Time mode, angle = " + (angle/(2*pi))*360 + " degrees");
    plot(C(:,1),C(:,2),"-","Color","black");
    scatter(C(:,1),C(:,2),".","filled","MarkerFaceColor","black","MarkerEdgeColor","black");

    if size(C,1) == length(labels)
        text(C(:,1),C(:,2),labels,"VerticalAlignment","bottom","HorizontalAlignment","right");
    elseif size(C,1) == length(labels2)
        text(C(:,1),C(:,2),labels2,"VerticalAlignment","bottom","HorizontalAlignment","right");
    else
        labels = NaN;
    end

elseif numFactors==3
    % Feature vectors
    subplot(1,3,1),hold on;
    xlabel("Feature vector 1");
    ylabel("Feature vector 2");
    zlabel("Feature vector 3");
    for i=1:max(feature_mode(:,4))
        plot3(feature_mode(feature_mode(:,4)==i,1),feature_mode(feature_mode(:,4)==i,2),feature_mode(feature_mode(:,4)==i,3),".", "MarkerSize", 20, "Color",FeatureColors(i,:));
    end
    view(45,45);

    % feature legend
    h = zeros(max(feature_mode(:,numFactors+1)),1);
    for feature=1:max(feature_mode(:,numFactors+1))
        h(feature) = plot(NaN, NaN, "s", "MarkerEdgeColor", FeatureColors(feature,:), "MarkerFaceColor", FeatureColors(feature,:));
    end
    legend(h, feature_names, "Location", 'southoutside', "Orientation", "horizontal", "NumColumns", 3);

    % Individual vectors
    subplot(1,3,2),hold on;
    xlabel("Individual vector 1");
    ylabel("Individual vector 2");
    zlabel("Individual vector 3");
    for i=0:2
        plot3(individual_mode(individual_mode(:,4) == i,1),individual_mode(individual_mode(:,4) == i,2),individual_mode(individual_mode(:,4) == i,3),".", "MarkerSize", 20, "Color",RFcolors(i+1,:));
    end
    view(45,45);

    % Add id legend
    h = zeros(3, 1);
    h(1) = plot(NaN,NaN,'s', "MarkerEdgeColor", RFcolors(1,:), "MarkerFaceColor", RFcolors(1,:));
    h(2) = plot(NaN,NaN,'s', "MarkerEdgeColor", RFcolors(2,:), "MarkerFaceColor", RFcolors(2,:));
    h(3) = plot(NaN,NaN,'s', "MarkerEdgeColor", RFcolors(3,:), "MarkerFaceColor", RFcolors(3,:));
    legend(h, {'Low response','Mid response','High response'}, "Location", 'southoutside', "Orientation", "horizontal");
    
    % Time vectors
    labels = string(1:7);
    subplot(1,3,3),hold on;
    xlabel("Time vector 1");
    ylabel("Time vector 2");
    zlabel("Time vector 3");
    plot3(C(:,1),C(:,2),C(:,3),"-*");
    text(C(:,1),C(:,2),C(:,3),labels,"VerticalAlignment","bottom","HorizontalAlignment","right");
    view(45,45);
end

sgtitle(plot_title);
saveas(gcf, path)
close()
end
