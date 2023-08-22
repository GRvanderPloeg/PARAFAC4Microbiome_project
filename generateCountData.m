% Generate count data
% Start from the very end, a finished PARAFAC model
addpath(".\Matlab scripts\Scripts\"); % own scripts
addpath(".\Matlab scripts\N-way toolbox\"); % from Rasmus Bro

%%
% Define properties of the simulation
A_relative_strength = 1.0;
B_relative_strength = 1.0;
C_relative_strength = 1.0;

% Subject loadings
Ameans = [1 -3; -7 10]; % rows is component, cols is RFgroup
Astds = [10 10; 10 10];   % rows is component, cols is RFgroup
Anum = [25 25];               % number of subjects per RFgroup

% Feature loadings
Bmeans = [3 -5 -2 1; 1 -1 -2 5]; % rows is component, cols is feature group
Bstds = [5 5 5 5; 5 5 5 5]; % rows is component, cols is feature group
Bnum = [25 25 25 25];

% Time loadings
C = [0.1 0.9; 0.11 0.8; 0.2 0.65; 0.55 0.55; 0.3 0.5; 0.25 0.44; 0.2 0.3; 0.11 0.25; 0.1 0.65];

% Residuals
residualNoisePerc = 0.0;

% Added stds
stdMean = 1.90;             % TIFN: 1.9
stdStd = 0.57;              % TIFN: 0.57

% Added means
meansMean = 4.00;           % TIFN: 4
meansStd = 1.72;            % TIFN: 1.72

% Added geoMeans
geoMeansMean = 1.0328;      % TIFN: 1.0328
geoMeansStd = 0.0082;       % TIFN: 0.0082

% Setting library size
percVariation = 0.2;
readsPerFeature = 100;

%%
% Simulate loadings
I = sum(Anum);
J = sum(Bnum);
K = size(C, 1);
numComponents = size(C, 2);

% Subject loadings and metadata
%A = normrnd(0, 11, [41 2]);
A = zeros(I, numComponents);
Ameta = zeros(I, numComponents);
Ameta(:,1) = 1:I; % subject "names"
iterator = 1;

for i=1:size(Anum, 2)
    numSubjects = Anum(i);
    for n=1:numSubjects
        A(iterator,:) = [normrnd(Ameans(1,i), Astds(1,i), 1, 1) normrnd(Ameans(2,i), Astds(2,i), 1, 1)];
        Ameta(iterator, 2) = i;
        iterator = iterator + 1;
    end
end

Ameta = string(Ameta);

% Feature loadings and metadata
%B = normrnd(0.075, 0.075, [78 2]);
B = zeros(J, numComponents);
Bmeta = zeros(J, numComponents);
Bmeta(:,1) = 1:J;
iterator = 1;

for i=1:size(Bnum, 2)
    numFeatures = Bnum(i);
    for n=1:numFeatures
        B(iterator,:) = [normrnd(Bmeans(1,i), Bstds(1,i), 1, 1) normrnd(Bmeans(2,i), Bstds(2,i), 1, 1)];
        Bmeta(iterator,2) = i;
        iterator = iterator + 1;
    end
end

Bmeta = string(Bmeta);

% Normalize to length 1
% for n=1:numComponents
%     A(:,n) = A(:,n) / norm(A(:,n));
%     B(:,n) = B(:,n) / norm(B(:,n));
%     C(:,n) = C(:,n) / norm(C(:,n));
% end

%%
% Create processed matrix from vectors
M = A * krb(C,B)';
Mcube = reshape(M, I, J, K);

% Add residuals
% error ~ N(0, set percentage of real value)
residuals = zeros(size(Mcube));

for i=1:I
    for j=1:J
        for k=1:K
            residuals(i,j,k) = normrnd(0, abs(Mcube(i,j,k))*residualNoisePerc, 1, 1);
        end
    end
end

Mcube = Mcube + residuals;

%%
% Reverse scaling
%fakeStds = microb_tongue_stds; %normrnd(1.9, 0.57, [78, 1]);
stds = normrnd(stdMean, stdStd, [J 1]);

Mcube_unScl = Mcube;

for j=1:J
    for k=1:K
        Mcube_unScl(:,j,k) = Mcube_unScl(:,j,k) * stds(j);
    end
end

%%
% Reverse centering
%fakeMeans = microb_tongue_means; %normrnd(4, 1.72, [78, 7]);
means = normrnd(meansMean, meansStd, [J K]);

Mcube_unScl_unCnt = Mcube_unScl;

for k=1:K
    for j=1:J
        Mcube_unScl_unCnt(:,j,k) = Mcube_unScl_unCnt(:,j,k) + means(j,k);
    end
end

%%
% Reverse CLR
%fakeGeoMeans = tongueGeoMeans;
geoMeans = normrnd(geoMeansMean, geoMeansStd, [I*K 1]);

dummy = permute(Mcube_unScl_unCnt, [1 3 2]);
Mfinal = reshape(dummy, I*K, J);

for i=1:(I*K)
    Mfinal(i,:) = exp(Mfinal(i,:)) * geoMeans(i);
end

Mfinal = round(Mfinal) - 1; % make into integers and remove pseudocount

%%
% Replace negative counts with zero
min(Mfinal, [], "all") % report minimum value as diagnostic
mask = Mfinal < 0;
Mfinal(mask) = 0;

%%
% Correct for absurd total counts per sample
libSizes = round(normrnd(J*readsPerFeature, J*readsPerFeature*percVariation, I*K, 1));

for i=1:(I*K)
    Mfinal(i,:) = (Mfinal(i,:) / sum(Mfinal(i,:))) * libSizes(i);
end

Mfinal = round(Mfinal);

%%
% Save data and metadata
writematrix(Mfinal, "./tempData.csv");
writematrix(Ameta, "subjectMetadata.csv");
writematrix(Bmeta, "featureMetadata.csv");
writematrix(A, "subjectLoadings.csv");
writematrix(B, "featureLoadings.csv");
writematrix(C, "timeLoadings.csv");
