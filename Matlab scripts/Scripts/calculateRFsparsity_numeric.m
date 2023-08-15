function [sparsityLow, sparsityMid, sparsityHigh]=calculateRFsparsity_numeric(X, rf_low, rf_mid, rf_high)

Xzeros = X(:, 6:end) == 0;
low_responders = ismember(X(:,2), rf_low);
mid_responders = ismember(X(:,2), rf_mid);
high_responders = ismember(X(:,2), rf_high);

[I,J] = size(Xzeros);
sparsityLow = zeros(J,1);
sparsityMid = zeros(J,1);
sparsityHigh = zeros(J,1);

for j=1:J
    dataVector = Xzeros(:,j);
    sparsityLow(j) = sum(dataVector(low_responders)) / sum(low_responders) * 100;
    sparsityMid(j) = sum(dataVector(mid_responders)) / sum(mid_responders) * 100;
    sparsityHigh(j) = sum(dataVector(high_responders)) / sum(high_responders) * 100;
end