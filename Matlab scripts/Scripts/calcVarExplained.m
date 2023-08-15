function varExplained=calcVarExplained(X, A, B, C)
M = A*krb(C,B)';
M = reshape(M, size(X));

varExplained = (sumsqr(M) / sumsqr(X)) * 100;