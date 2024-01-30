function [Xhat, PRESS] = crossValPCA(X,maxF)
[I,J] = size(X);

% Initialize output matrices
Xhat = cell(1,maxF);
for f=1:maxF
    Xhat{f} = zeros(size(X));
end

% For each sample i
for i=1:I
    % Remove sample i
    X_rm = X;
    X_rm(i,:) = [];

    % Calculate PCA on the remainder
    [U, S, V] = svd(X_rm, "econ");
    T = zeros(size(U));
    variance = diag(S);
    for n=1:size(S,2)
        T(:,n) = U(:,n) * variance(n,:);
    end
    P = V;

    % For each factor
    for f=1:maxF
        Tf = T(:,f);
        Pf = P(:,f);

        % For every variable
        for j=1:J
            P_minusJ = Pf;
            P_minusJ(j) = [];
            
            x_i_minusJ = X(i,:);
            x_i_minusJ(j) = [];
    
            % Estimate the score
            t_minusJ = x_i_minusJ * P_minusJ * inv(P_minusJ' * P_minusJ);
    
            % Estimate the element
            Xhat{f}(i,j) = t_minusJ * Pf(j);
        end
    end
end

PRESS = zeros(1,f);
for f=1:maxF
    error = X;
    for n=1:f
        error = error - Xhat{n};
    end
    PRESS(f) = sumsqr(error);
end
