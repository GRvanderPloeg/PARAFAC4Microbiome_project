function plotCORCONDIA(X, maxFactors, numReps, path)

Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 2;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)

core_consistency_mean = 1:maxFactors;
core_consistency_std = 1:maxFactors;

for i=1:maxFactors
    dummy = 1:numReps;
    for j=1:numReps
        [~, ~, ~, c] = myParafac(X, i, 50, Options);
        dummy(j) = c;
    end

    core_consistency_mean(i) = mean(dummy);
    core_consistency_std(i) = std(dummy);
end

subplot(1,1,1), errorbar(1:maxFactors, core_consistency_mean, core_consistency_std)
ylim([0, 100])
xlabel("Number of components")
ylabel("Core consistency")

saveas(gcf, path)

end