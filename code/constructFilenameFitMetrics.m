function [filenameFitMetricsOutput] = constructFilenameFitMetrics(...
    isMeansOfMeans,isLogTransformed,isFluxNormalised,choiceZref)

    baseName = 'fitmetrics_classic_';

    % Parts for filename 
    infix1 = {'annualfits','meansofmeans'};
    infix2 = {'notlog','log'};
    infix3 = {'notnormalised','normalised'};
    infix4 = {'100m','zeu','inflexion'};

    % Use logical indexing to set infix values
    infix1 = infix1{isMeansOfMeans + 1};
    infix2 = infix2{isLogTransformed + 1};
    infix3 = infix3{isFluxNormalised + 1};
    infix4 = infix4{choiceZref};

    % Combine all parts into the filename
    filenameFitMetricsOutput = sprintf('%s%s_%s_%s_%s.mat',...
        baseName,infix1,infix2,infix3,infix4);

end