function [Xcube,universalSubjects]=rawDataToCube(Xlong, Xmeta_subjects, Xmeta_timepoints, numTimepoints)

subjectNames = unique(Xmeta_subjects);
universalSubjects = false(length(subjectNames), 1);
for i=1:length(subjectNames)
    name = subjectNames(i);
    numberOfSamples = sum(Xmeta_subjects == name);
    universalSubjects(i) = (numberOfSamples == numTimepoints);
end
universalSubjects = subjectNames(universalSubjects);

Xmeta = [Xmeta_subjects Xmeta_timepoints];
Xmeta(:,end+1) = 1:size(Xmeta,1);
Xmeta = sortrows(Xmeta, [2 1]); % sort on visit number (ascending), then subject name (alphabetical)
Xmeta = Xmeta(ismember(Xmeta(:,1), universalSubjects), :);
keepRowIndices = str2double(Xmeta(:,end));

Xlong = Xlong(keepRowIndices, :);
I = size(Xlong, 1) / numTimepoints;
J = size(Xlong, 2);
K = numTimepoints;

Xcube = reshape(Xlong, I, K, J);
Xcube = permute(Xcube, [1 3 2]);

