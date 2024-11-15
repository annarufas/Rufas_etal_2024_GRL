function [selectedDepths,selectedFluxes,selectedErrs] = extractDataBelowZref(...
    depthVector,fluxVector,fluxErrVector,choiceZref,locZeu,maxZref)

% EXTRACTDATABELOWZREF Extracts POC flux values and associated depths below
% the reference depth (zref).
%
%   INPUT: 
%       depthVector   - vector of depth values
%       fluxVector    - vector of POC flux values
%       fluxErrVector - vector of POC flux errors
%       choiceZref    - choice of reference depth (zeu, 100 m or inflexion point)
%       locZeu        - local euphotic layer depth
%       maxZref       - maximum allowable reference depth
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

validIdx = ~isnan(depthVector); 
validDepths = depthVector(validIdx);
validFluxes = fluxVector(validIdx);
validErrs = fluxErrVector(validIdx);

%% Initialise outputs

selectedDepths = [];
selectedFluxes = [];
selectedErrs = [];

%% Determine zref based on choiceZref

switch choiceZref
    case 1  % Fixed 100 m depth
        [~, idxTargetDepth] = min(abs(validDepths - 100)); 
        
    case 2  % Local euphotic depth (zeu)
        [~, idxTargetDepth] = min(abs(validDepths - locZeu)); 
        
    case 3  % Depth at inflection point 
        % Find the depth where flux first decreases (simple check for inflection)
        idxTargetDepth = [];
        for i = 1:length(validFluxes)-1
            if all(validFluxes(i) > validFluxes(i+1:end))
                idxTargetDepth = i;
                break;
            end
        end

end

%% Extract

if ~isempty(idxTargetDepth)
    % Filter depths greater than or equal to zref (target depth)
    selectedDepths = validDepths(validDepths >= validDepths(idxTargetDepth));
    
    % Ensure that selected depths do not exceed maxZref
    if selectedDepths(1) <= maxZref
        selectedFluxes = validFluxes(validDepths >= validDepths(idxTargetDepth)); 
        selectedErrs = validErrs(validDepths >= validDepths(idxTargetDepth)); 
    else
        selectedDepths = [];
        selectedFluxes = [];
        selectedErrs = [];
    end
end

end % extractDataBelowZref