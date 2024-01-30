function numComponentsOverviewPlot(X, maxFactors, numReps, numSegments, overallTitle, path)
I = size(X,1); % Number of rows
J = size(X,2); % Number of columns
K = size(X,3); % Number of timepoints

Options = [];
Options(1) = 1e-3;  % convergence (default 1e-6)
Options(2) = 2;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)
Options(4) = 2;     % no scaling applied (default 1)
Options(5) = NaN;   % no output given (default 10)
Options(6) = 1000;  % max number of iterations (default 2500)

CV_options = [];
CV_options(1) = 1e-3;
CV_options(2) = 1;
CV_options(3) = 0;
CV_options(4) = 2;     % no scaling applied (default 1)
CV_options(5) = NaN;   % no output given (default 10)
CV_options(6) = 1000;  % max number of iterations (default 2500)

resSSQ = zeros(numReps,maxFactors);
corcondia = zeros(numReps,maxFactors);
numIterations = zeros(numReps,maxFactors);
PRESS = zeros(maxFactors,1);
tuckerCongruence = zeros(numReps,maxFactors);
varExp_full = zeros(numReps,maxFactors);
varExp_CV = zeros(maxFactors,1);

for i=1:maxFactors
    disp("Working on " + i + "-component models.")

    for j=1:numReps
        % Fit a PARAFAC model and note the statistics
        [model, iterations, error, corcon, ~,~,~] = mySilentParafac(X,i,1,Options);
        [A, B, C] = fac2let(model);
        resSSQ(j,i) = error;
        corcondia(j,i) = corcon;
        numIterations(j,i) = iterations;
        varExp_full(j,i) = 100 - (error / sumsqr(X)) * 100;
        tuckerCongruence(j,i) = calcTuckerCongruence(A,B,C);
    end
    
    % Crossvalidation of PRESS + variance explained
    CVresiduals = zeros(I, J*K);
    for k=1:numSegments
        Xpart = reshape(X, I, J*K);
        Xpart(k:numSegments:end) = NaN;
        Xpart = reshape(Xpart, I, J, K);
    
        CVmodel = mySilentParafac(Xpart,i,1,CV_options);
        [A,B,C] = fac2let(CVmodel);
        modeled_data = A*krb(C,B)';
        modeled_data = reshape(modeled_data, I, J, K);
        residuals = reshape((X - modeled_data), I, J*K);
        CVresiduals(k:numSegments:end) = residuals(k:numSegments:end);
    end
        PRESS(i) = sumsqr(CVresiduals);
        varExp_CV(i) = 100 - (sumsqr(CVresiduals) / sumsqr(X)) * 100;
end

set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])

subplot(3,2,1); errorbar(1:maxFactors, mean(resSSQ), std(resSSQ)); hold on; plot(1:maxFactors, PRESS); yline(sumsqr(X),'r--', "Input data SSQ", 'LineWidth', 2); title("Sum of squares"); ylim([0 sumsqr(X)*1.2]);
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'s', "MarkerEdgeColor", [0, 0.4470, 0.7410], "MarkerFaceColor", [0, 0.4470, 0.7410]);
h(2) = plot(NaN,NaN,'s', "MarkerEdgeColor", [0.8500, 0.3250, 0.0980], "MarkerFaceColor", [0.8500, 0.3250, 0.0980]);
legend(h, {'Residual','PRESS (CV)'}, "Location", 'southoutside', "Orientation", "horizontal");

subplot(3,2,2); errorbar(1:maxFactors, mean(corcondia), std(corcondia)); title("Core consistency"); ylim([0 100]);
subplot(3,2,3); errorbar(1:maxFactors, mean(numIterations), std(numIterations)); title("Number of iterations");
subplot(3,2,4); plot(1:maxFactors, mean(varExp_full)); hold on; plot(1:maxFactors, varExp_CV); title("Variance explained"); ylim([0 100]);
subplot(3,2,5); errorbar(2:maxFactors, mean(tuckerCongruence(:,2:end)), std(tuckerCongruence(:,2:end))); title("Tucker congruence"); ylim([0 1]);

rawSSQ = zeros(1,maxFactors) + sumsqr(X);
subplot(3,2,6); bar([rawSSQ; mean(resSSQ); (PRESS')]'); hold on; title("SSQ bar plot"); ylim([0 sumsqr(X)*1.2]);
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'s', "MarkerEdgeColor", [0, 0.4470, 0.7410], "MarkerFaceColor", [0, 0.4470, 0.7410]);
h(2) = plot(NaN,NaN,'s', "MarkerEdgeColor", [0.8500, 0.3250, 0.0980], "MarkerFaceColor", [0.8500, 0.3250, 0.0980]);
h(3) = plot(NaN,NaN,'s', "MarkerEdgeColor", [0.9290, 0.6940, 0.1250], "MarkerFaceColor", [0.9290, 0.6940, 0.1250]);
legend(h, {'Input data','Residual','PRESS (CV)'}, "Location", 'southoutside', "Orientation", "horizontal");


sgtitle(overallTitle);
saveas(gcf, path)
close()
end