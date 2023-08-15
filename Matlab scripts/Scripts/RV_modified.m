function [RV]=RV_modified(X,Y)
AA=X*X';
BB=Y*Y';
AA0 = AA - diag(diag(AA),0);
BB0 = BB - diag(diag(BB),0);
RV = trace(AA0*BB0)/sumsqr(AA0)^.5/sumsqr(BB0)^.5;

