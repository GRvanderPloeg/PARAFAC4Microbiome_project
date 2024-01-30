function [coefficient_means, coefficient_stds]=CV_coeff_NPLS(X, Y, numComponents, mode)

    numSubjects = size(X, 1);
    numFeatures = size(X, 2);
    percCV = 10; % percentage of samples to cut out if not using LOO

    if mode==1
        % Leave one out CV
        result = cell(numComponents, 1);
    
        for i=1:numSubjects
            X_train = X;
            X_train(i,:,:) = []; % remove one subject
            Y_train = Y;
            Y_train(i) = []; % remove one subject
    
            [Xfactors,~,~,~,~,~,~,reg] = npls(X_train, Y_train, numComponents);
            coefficients = calculateCoefficientsNPLS(Xfactors, reg);
    
            for j=1:numComponents
                result{j}(:,i) = coefficients(:,j);
            end
        end
    elseif mode == 2
        % X percent CV
        cutSamples = fix(numSubjects/percCV);
        numReps = fix(numSubjects/cutSamples);
        result = cell(numComponents, numReps);
    
        for i=1:numReps
            cutOut = ((4*i)-3):(4*i);
            X_train = X;
            X_train(cutOut,:,:) = []; % remove one subject
            Y_train = Y;
            Y_train(cutOut) = []; % remove one subject
    
            [Xfactors,~,~,~,~,~,~,reg] = npls(X_train, Y_train, numComponents);
            coefficients = calculateCoefficientsNPLS(Xfactors, reg);
    
            for j=1:numComponents
                result{j}(:,i) = coefficients(:,j);
            end
        end
    end

    coefficient_means = zeros(numFeatures, numComponents);
    coefficient_stds = zeros(numFeatures, numComponents);
    
    for i=1:numComponents
        coefficient_means(:,i) = mean(result{i}, 2);
        coefficient_stds(:,i) = std(result{i}, 0, 2);
    end
end
