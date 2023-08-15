function savePARAFAC(X, model, individual_mode_metadata, feature_mode_metadata, path)

[A, B, C] = fac2let(model);
[A, B, C] = sortParafacComponents(X, A, B, C);
num_components = size(A,2);

individual_mode = [A individual_mode_metadata];
feature_mode = [B feature_mode_metadata];
time_mode = C;

M = A*krb(C,B)';
M = reshape(M, size(X));

writematrix(individual_mode, path + "_individual_mode.csv");
writematrix(feature_mode, path + "_feature_mode.csv");
writematrix(time_mode, path + "_time_mode.csv");
writematrix(M, path + "_model.csv");
writematrix(X, path + "_input.csv");

for i=1:num_components
    Mcomp = A(:,i)*krb(C(:,i),B(:,i))';
    Mcomp = reshape(Mcomp, size(X));
    writematrix(Mcomp, path + "_component_" + i + ".csv")
end