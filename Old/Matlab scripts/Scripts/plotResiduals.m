function plotResiduals(X, numComponents, numTimepoints, path)
    [Factors] = parafac(X,numComponents);
    M = nmodel(Factors);
    E = X - M;
    
    plotDim = ceil(sqrt(numTimepoints));
    for i=1:numTimepoints
        subplot(plotDim,plotDim,i),mesh(X(:,:,i))
    end
    saveas(gcf, path+"_raw.jpg")
    
    for i=1:numTimepoints
        subplot(plotDim,plotDim,i),mesh(E(:,:,i))
    end
    saveas(gcf, path+"_residuals.jpg")
end