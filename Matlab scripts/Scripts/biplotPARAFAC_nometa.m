function biplotPARAFAC_nometa(X, model, featureOffset, plot_title, path)
% Only supports 2 and 3-component PARAFAC models!

% FeatureOffset = 1 for metabolomics (if you want super_pathway colors)
% = 2 for microbiome (if you want phylum colors)

if ~exist('path','var')
    path="";
end

[A, B, C] = fac2let(model);
[A, B, C] = sortParafacComponents(X, A, B, C);
numFactors = size(A,2);

individual_mode = A;
feature_mode = B;

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
    plot(Vector1, Vector2, ".", "MarkerSize", 20);

    % Individual vectors
    subplot(1,3,2),hold on;
    xlabel("Individual vector 1");
    ylabel("Individual vector 2");
    Vector1 = individual_mode(:,1);
    Vector2 = individual_mode(:,2);
    angle = acos((Vector1'*Vector2)/(norm(Vector1)*norm(Vector2)));
    title("Individual mode, angle = " + (angle/(2*pi))*360 + " degrees");
    plot(individual_mode(:,1), individual_mode(:,2),".", "MarkerSize", 20);
    
    % Time vectors
    subplot(1,3,3),hold on;
    xlabel("Time vector 1");
    ylabel("Time vector 2");
    Vector1 = C(:,1);
    Vector2 = C(:,2);
    angle = acos((Vector1'*Vector2)/(norm(Vector1)*norm(Vector2)));
    title("Time mode, angle = " + (angle/(2*pi))*360 + " degrees");
    plot(C(:,1),C(:,2),"-","Color","black");
    scatter(C(:,1),C(:,2),".","filled","MarkerFaceColor","black","MarkerEdgeColor","black");
    labels=string(1:size(C,1));
    text(C(:,1),C(:,2),labels,"VerticalAlignment","bottom","HorizontalAlignment","right");
end

sgtitle(plot_title);

if path ~= ""
    saveas(gcf, path)
    close()
end

