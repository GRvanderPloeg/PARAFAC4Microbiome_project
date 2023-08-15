function benchmarkParafac(X, model, realI, realJ, realK, overallTitle, path)

if ~exist('path','var')
    path="";
end

I = size(X,1);
J = size(X,2);
K = size(X,3);
colormap = lines();

[A,B,C] = fac2let(model);
numComponents_model = size(A,2);
numComponents_real = size(realI,2);

% Remove variance from A to make loadings comparable
A_novar = zeros(size(A));
for n=1:numComponents_model
    A_novar(:,n) = A(:,n) ./ sqrt(sum(A(:,n).^2));
end

if numComponents_model == numComponents_real
    [A_plot, B_plot, C_plot] = mapRealModeledVectors(realI, realJ, realK, A_novar, B, C);
else
    A_plot = A_novar;
    B_plot = B;
    C_plot = C;
end

M = A*krb(C,B)';
M = reshape(M, size(X));

set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])

% Plot 1: true loadings of the 'individuals'
subplot(3,3,1); bar(realI); title("Real 'individual' loadings"); hold on;

% Legend
h = zeros(numComponents_real, 1);
labels = cell(1,numComponents_real);
for n=1:numComponents_real
    col = colormap(n,:);
    h(n) = plot(NaN,NaN,'s', "MarkerEdgeColor", col, "MarkerFaceColor", col);
    labels{n} = "True loading " + n;
end
legend(h, labels, "Location", 'southoutside', "Orientation", "horizontal");

% Plot 2 modeled loadings of the 'individuals'
subplot(3,3,2); bar(A_plot); title("Modeled 'individual' loadings"); hold on;

% Legend
h = zeros(numComponents_model, 1);
labels = cell(1,numComponents_model);
for n=1:numComponents_model
    col = colormap(n,:);
    h(n) = plot(NaN,NaN,'s', "MarkerEdgeColor", col, "MarkerFaceColor", col);
    labels{n} = "Modeled loading " + n;
end
legend(h, labels, "Location", 'southoutside', "Orientation", "horizontal");

% Plot 3: ssq of the 'individuals'
ssq_individuals = zeros(I, 2);
for i=1:I
    ssq_individuals(i,1) = sumsqr(X(i,:,:));
    ssq_individuals(i,2) = sumsqr(M(i,:,:));
end
subplot(3,3,3); bar(ssq_individuals); title("SSQ per 'individual'"); hold on;
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'s', "MarkerEdgeColor", [0, 0.4470, 0.7410], "MarkerFaceColor", [0, 0.4470, 0.7410]);
h(2) = plot(NaN,NaN,'s', "MarkerEdgeColor", [0.8500, 0.3250, 0.0980], "MarkerFaceColor", [0.8500, 0.3250, 0.0980]);
legend(h, {'Input data','Model'}, "Location", 'southoutside', "Orientation", "horizontal");

% Plot 4: true loadings of the 'features'
subplot(3,3,4); bar(realJ); title("Real 'feature' loadings"); hold on;

% Legend
h = zeros(numComponents_real, 1);
labels = cell(1,numComponents_real);
for n=1:numComponents_real
    col = colormap(n,:);
    h(n) = plot(NaN,NaN,'s', "MarkerEdgeColor", col, "MarkerFaceColor", col);
    labels{n} = "True loading " + n;
end
legend(h, labels, "Location", 'southoutside', "Orientation", "horizontal");

% Plot 5: modeled loadings of the 'features'
subplot(3,3,5); bar(B_plot); title("Modeled 'feature' loadings"); hold on;

% Legend
h = zeros(numComponents_model, 1);
labels = cell(1,numComponents_model);
for n=1:numComponents_model
    col = colormap(n,:);
    h(n) = plot(NaN,NaN,'s', "MarkerEdgeColor", col, "MarkerFaceColor", col);
    labels{n} = "Modeled loading " + n;
end
legend(h, labels, "Location", 'southoutside', "Orientation", "horizontal");

% Plot 6: ssq of the 'features'
ssq_features = zeros(J, 2);
for j=1:J
    ssq_features(j,1) = sumsqr(X(:,j,:));
    ssq_features(j,2) = sumsqr(M(:,j,:));
end
subplot(3,3,6); bar(ssq_features); title("SSQ per 'feature'"); hold on;
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'s', "MarkerEdgeColor", [0, 0.4470, 0.7410], "MarkerFaceColor", [0, 0.4470, 0.7410]);
h(2) = plot(NaN,NaN,'s', "MarkerEdgeColor", [0.8500, 0.3250, 0.0980], "MarkerFaceColor", [0.8500, 0.3250, 0.0980]);
legend(h, {'Input data','Model'}, "Location", 'southoutside', "Orientation", "horizontal");

% Plot 7: real loadings of 'time'
subplot(3,3,7); plot(1:K, realK); title("Real 'time' mode"); hold on;

% Legend
h = zeros(numComponents_real, 1);
labels = cell(1,numComponents_real);
for n=1:numComponents_real
    col = colormap(n,:);
    h(n) = plot(NaN,NaN,'s', "MarkerEdgeColor", col, "MarkerFaceColor", col);
    labels{n} = "True loading " + n;
end
legend(h, labels, "Location", 'southoutside', "Orientation", "horizontal");

% Plot 8: modeled loadings of 'time'
subplot(3,3,8); plot(1:K, C_plot); title("Modeled 'time' mode"); hold on;

% Legend
h = zeros(numComponents_model, 1);
labels = cell(1,numComponents_model);
for n=1:numComponents_model
    col = colormap(n,:);
    h(n) = plot(NaN,NaN,'s', "MarkerEdgeColor", col, "MarkerFaceColor", col);
    labels{n} = "Modeled loading " + n;
end
legend(h, labels, "Location", 'southoutside', "Orientation", "horizontal");

% Plot 9: overall ssq
xvalues = categorical({'Input data', 'Model'});
xvalues = reordercats(xvalues,{'Input data', 'Model'});
subplot(3,3,9); title("SSQ overall"); hold on;
b = bar(xvalues, [sumsqr(X), sumsqr(M)]);
b.FaceColor = 'flat';
b.CData(2,:) = [0.8500, 0.3250, 0.0980];

h = zeros(2, 1);
h(1) = plot(NaN,NaN,'s', "MarkerEdgeColor", [0, 0.4470, 0.7410], "MarkerFaceColor", [0, 0.4470, 0.7410]);
h(2) = plot(NaN,NaN,'s', "MarkerEdgeColor", [0.8500, 0.3250, 0.0980], "MarkerFaceColor", [0.8500, 0.3250, 0.0980]);
legend(h, {'Input data','Model'}, "Location", 'southoutside', "Orientation", "horizontal");

sgtitle(overallTitle);

if path ~= ""
    saveas(gcf, path)
    close()
end
