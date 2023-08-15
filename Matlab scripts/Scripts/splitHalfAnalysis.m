function splitHalfAnalysis(X, numComponents, path)
    microb_tongue_splitA = X(1:round(size(X,1)/2),:,:);
    microb_tongue_splitB = X(round(size(X,1)/2)+1:size(X,1),:,:);
    microb_tongue_splitC = X(1:2:end,:,:);
    microb_tongue_splitD = X(2:2:end,:,:);

    randomDraw = randperm(size(X,1));
    splitE = randomDraw(1:round(size(X,1)/2));
    splitF = randomDraw(round(size(X,1)/2)+1:size(X,1));
    microb_tongue_splitE = X(splitE,:,:);
    microb_tongue_splitF = X(splitF,:,:);
    
    plotPARAFAC(microb_tongue_splitA, numComponents, 250, path+"_splitA.jpg");
    plotPARAFAC(microb_tongue_splitB, numComponents, 250, path+"_splitB.jpg");
    plotPARAFAC(microb_tongue_splitC, numComponents, 250, path+"_splitC.jpg");
    plotPARAFAC(microb_tongue_splitD, numComponents, 250, path+"_splitD.jpg");
    plotPARAFAC(microb_tongue_splitE, numComponents, 250, path+"_splitE.jpg");
    plotPARAFAC(microb_tongue_splitF, numComponents, 250, path+"_splitF.jpg");
end