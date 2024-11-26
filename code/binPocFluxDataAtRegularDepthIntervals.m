function [depthsBinned,fluxBinned] = binPocFluxDataAtRegularDepthIntervals(...
    profileDepths,profileAvg,profileErrTot,nData,maxNumDepths)

% BINPOCFLUXDATAATREGULARDEPTHINTERVALS Bins particulate organic carbon (POC) 
% flux data at regular depth intervals and propagates associated errors. 
% This binning optimises the data for smoother curve fitting by averaging 
% flux values across consistent depth bins.
%
%   INPUT: 
%       profileDepths - vector of depth values
%       profileAvg    - vector of POC flux values
%       profileErrTot - vector of POC flux errors
%       nData         - count of samples used to derive each mean POC flux value in profileAvg
%       maxNumDepths  - parameter
%
%   OUTPUT:
%       depthsBinned - vector of depth centres for each bin
%       fluxBinned   - vector of binned POC flux values and their errors
%
%   This script uses this external function: 
%       worstcase.m  - FileExchange
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 13 Nov 2024  
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------
                
%% Set bin size and edges

binSize = 5;
maxDepth = 5000;

%% Allocate memory for binned data

depthsBinned = NaN(maxNumDepths,1);
fluxBinned = NaN(maxNumDepths,2); % 2nd dim: 1=avg, 2=err

%% Build the bins

firstDepth = profileDepths(1);
roundedDepth = round(firstDepth / binSize) * binSize;
binEdges = (roundedDepth:binSize:maxDepth)';
[~, ~, idxBin] = histcounts(profileDepths, binEdges);

%% Bin data

fluxAvgBinned = NaN(length(binEdges), 1);
fluxErrBinned = NaN(length(binEdges), 1);
depthCentre = NaN(length(binEdges), 1);

for iBin = 1:length(binEdges)

    % Values within the current bin
    inBin = (idxBin == iBin);
    fluxBinMembers = profileAvg(inBin);
    errBinMembers = profileErrTot(inBin);
    nSamples = nData(inBin);

    if sum(inBin) > 1
        
        % Weighted average calculation
        paramsWeightedAverage = fluxBinMembers .* nSamples;
        calculateWeightedAvg = @(x) sum(x) ./ sum(nSamples(:), 'omitnan');
        fluxAvgBinned(iBin) = calculateWeightedAvg(paramsWeightedAverage);

        % Calculate propagated error (worst-case)
        [~,~,~,f_MID,f_UB,~,~] = worstcase(...
            @(paramsWeightedAverage) calculateWeightedAvg(paramsWeightedAverage),...
            paramsWeightedAverage, errBinMembers);
        fluxErrBinned(iBin) = f_UB - f_MID;
        depthCentre(iBin) = binEdges(iBin);
    
    % Nothing to average
    elseif sum(inBin) == 1
        
        fluxAvgBinned(iBin) = fluxBinMembers;
        fluxErrBinned(iBin) = errBinMembers;
        depthCentre(iBin) = binEdges(iBin);
        
    end
end

% Filter out NaNs and assign binned values
validIdx = ~isnan(fluxAvgBinned);
nValid = sum(validIdx);
fluxBinned(1:nValid,1) = fluxAvgBinned(validIdx);
fluxBinned(1:nValid,2) = fluxErrBinned(validIdx);
depthsBinned(1:nValid) = depthCentre(validIdx);
    
end