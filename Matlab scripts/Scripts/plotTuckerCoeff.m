function plotTuckerCoeff(X, maxFactors, path)

Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 2;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)

% Make matrix splittable by making the number of samples even
% If there is an uneven number, one sample randomly gets removed
if (mod(size(X,1),2) ~= 0)
    removedSample = randperm(size(X,1),1);
    selectedSamples = [1:removedSample-1 removedSample+1:size(X,1)];
    X = X(selectedSamples,:,:);

% Ordinary splitting in half
X_splitA = X(1:round(size(X,1)/2),:,:);
X_splitB = X(round(size(X,1)/2)+1:size(X,1),:,:);

tucker_coeff = 1:maxFactors;

for i=1:maxFactors
    [FactorsA,~,~,~,~,~,~,~] = myParafac(X_splitA, i, 250, Options);
    [FactorsB,~,~,~,~,~,~,~] = myParafac(X_splitB, i, 250, Options);

    [Multiphi,~] = ncosine(FactorsA,FactorsB);
    Multiphi
    tucker_coeff(i) = mean(max(Multiphi,[],2));
end

subplot(1,1,1), plot(1:maxFactors, tucker_coeff)
xlabel("Number of components")
ylabel("Tucker coefficient")

saveas(gcf, path)
end