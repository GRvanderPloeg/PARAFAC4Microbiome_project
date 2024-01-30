function result=logisticTransformation(X)
result = exp(X) ./ (1 + exp(X));
result = result - 0.5;
result = result ./ max(max(max(result)));