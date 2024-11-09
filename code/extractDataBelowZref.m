function [selectedDepths,selectedFluxes,selectedErrs] = extractDataBelowZref(...
    depthVector,fluxVector,fluxErrVector,choiceZref,locZeu)

% EXTRACTDATABELOWZREF Extracts POC flux values and associated depths below
% the reference depth (zref).
%
%   INPUT: 
%       depthVector   - vector of depth values
%       fluxVector    - vector of POC flux values
%       fluxErrVector - vector of POC flux errors
%       choiceZref    - choice of reference depth (zeu, 100 m or inflexion point)
%       locZeu        - local euphotic layer depth
%
%   OUTPUT:
%       selectedDepths - vector of selected depths
%       selectedFluxes - vector of selected POC flux values
%       selectedErrs   - vector of selected POC flux errors
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 4 Nov 2024  
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Get rid of NaNs

idxsNonNanDepths = find(~isnan(depthVector));
nonNanDepths = depthVector(idxsNonNanDepths);
nonNanFluxes = fluxVector(idxsNonNanDepths);
nonNanErrs = fluxErrVector(idxsNonNanDepths);

%% Initialise outputs

selectedDepths = [];
selectedFluxes = [];
selectedErrs = [];

%% Determine zref based on choiceZref

if (choiceZref == 1)
    [~,idxTargetDepth] = min(abs(nonNanDepths - 100));
elseif (choiceZref == 2)  
    [~,idxTargetDepth] = min(abs(nonNanDepths - locZeu));
elseif (choiceZref == 3)
    % Loop through each element to find the first point where 
    % all subsequent values are smaller
    idxTargetDepth = [];
    for i = 1:length(nonNanFluxes)-1
        if all(nonNanFluxes(i) > nonNanFluxes(i+1:end))
            idxTargetDepth = i;
            break;
        end
    end
end

%% Extract

if ~isempty(idxTargetDepth)
    idxAccDepths = nonNanDepths >= nonNanDepths(idxTargetDepth);
    selectedDepths = nonNanDepths(idxAccDepths);
    selectedFluxes = nonNanFluxes(idxAccDepths); 
    selectedErrs = nonNanErrs(idxAccDepths); 
end
    
end % extractDataBelowZref