function plotPARAFAC(X, numFactors, numReps, path)

Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 1;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)

[Factors, ~, ~, ~] = myParafac(X, numFactors, numReps, Options);
[A, B, C] = fac2let(Factors);

for i=0:(numFactors-1)
    varExplained = calcVarExplained(X, A(:,i+1), B(:,i+1), C(:,i+1));
    subplot(numFactors, 3, 1+i*3), bar(B(:,i+1)), ylabel("Comp. " + string(i+1) + " (" + string(varExplained) + "%)")
    subplot(numFactors, 3, 2+i*3), bar(A(:,i+1))
    subplot(numFactors, 3, 3+i*3), plot(C(:,i+1)), xlim([1, 7])
end

saveas(gcf, path)
end