function [P,Q,R]=pseudoRank(X)
% pseudoRank.m calculates the pseudo rank of three-way matrices.
% NOTE: DO NOT SUPPLY CENTERED/SCALED DATA.
%
% This algorithm computes every element x_ij per factor f of PCA to
% separate the noise from the variation. Make sure that you have done no
% operations on the data that connects x_hat_ij to x_ij. They have to be
% completely independent!
%
% Inputs:
%   X   : Raw 3-way input matrix of data
%
% Outputs:
%   P   : pseudo rank of mode I
%   Q   : pseudo rank of mode J
%   R   : pseudo rank of mode K

[I,J,K] = size(X);

% Calculate maximally allowed number of factors per direction
maxF_I = min(I-2, 500);
maxF_J = min(J-1, 500);
maxF_K = min(K-1, 500);

% Mode I
Xwide = reshape(X, I, J*K);
[~, PRESS_I] = crossValPCA(Xwide, maxF_I);
P = find(PRESS_I == min(PRESS_I));

% Mode J
Xperm = permute(X, [2 1 3]);
Xlong = reshape(Xperm, J, I*K);
[~, PRESS_J] = crossValPCA(Xlong, maxF_J);
Q = find(PRESS_J == min(PRESS_J));

% Mode 
Xperm = permute(X, [3 1 2]);
Xdeep = reshape(Xperm, K, I*J);
[~, PRESS_K] = crossValPCA(Xdeep, maxF_K);
R = find(PRESS_K == min(PRESS_K));

subplot(1,3,1); plot(PRESS_I);
subplot(1,3,2); plot(PRESS_J);
subplot(1,3,3); plot(PRESS_K);