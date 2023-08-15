function [mappedA, mappedB, mappedC]=mapRealModeledVectors(realI, realJ, realK, A, B, C)
numComponents = size(A,2);

mappedA = zeros(size(A,1), size(A,2));
mappedB = zeros(size(B,1), size(B,2));
mappedC = zeros(size(C,1), size(C,2));

modeledVectors = cell(1,3);
modeledVectors{1} = A;
modeledVectors{2} = B;
modeledVectors{3} = C;

realVectors = cell(1,3);
realVectors{1} = realI;
realVectors{2} = realJ;
realVectors{3} = realK;

result = cell(1,3);
result{1} = mappedA;
result{2} = mappedB;
result{3} = mappedC;

for x=1:3
    modeledVector = modeledVectors{x};
    realVector = realVectors{x};
    solution = result{x};
    disallowed = zeros(numComponents,1); % keep track of mapped vectors

    for n=1:numComponents
        v_model = modeledVector(:,n);
        bestError = 1e20;
        mapping = 0;
        inverse = 1;
    
        for m=1:numComponents
            v_real = realVector(:,m);
    
            if sqrt(sumsqr(v_real-v_model)) < bestError && (~ismember(m, disallowed))
                bestError = sumsqr(v_real-v_model);
                mapping = m;
                inverse = 1;
            end
    
            if sqrt(sumsqr(v_real-(-1*v_model))) < bestError && (~ismember(m, disallowed))
                bestError = sumsqr(v_real-(-1*v_model));
                mapping = m;
                inverse = -1;
            end
        end

        realVector(:,mapping)
        modeledVector(:,n)
        solution(:,mapping) = inverse .* modeledVector(:,n);
        disallowed(n) = mapping;

    end

    result{x} = solution;
  
end

mappedA = result{1};
mappedB = result{2};
mappedC = result{3};
