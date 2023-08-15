function coefficients=calculateCoefficientsNPLS(Xfactors, reg)
    [~, ~, C] = fac2let(Xfactors);
    numComponents = size(C,2);
    numCoefficients = size(reg{1}, 1);
    coefficients = zeros(numCoefficients, numComponents);

    for i=1:size(C,2)
        coefficients_comp = reg{i}' ./ C(:,i); 
        coefficients_comp = coefficients_comp(1,:)';
        coefficients(:,i) = coefficients_comp;
    end
end

    
