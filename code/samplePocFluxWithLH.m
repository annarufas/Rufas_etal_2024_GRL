function [paramMatrixFlux,paramMatrixDepth] = samplePocFluxWithLH(profileFlux_mean,...
    profileFlux_error,profileDepths,nSamples)

% SAMPLEPOCFLUXWITHLH Latin Hypercube sampling on POC flux data.
% This function generates a matrix of sampled POC flux values based on the
% mean flux values and their associated errors at specified depth points. 
% The resulting samples will be used for calculating Martin's b and the 
% remineralisation depth scale (z*).
% 
%   INPUT: 
%       profileFlux_mean  - POC flux cast mean values (mg C m-2 d-1)
%       profileFlux_error - POC flux cast errors (mg C m-2 d-1)
%       profileDepths     - POC flux cast depths (m)
%       nSamples          - number of random profiles we want to generate
%
%   OUTPUT:
%       paramMatrixFlux   - nFluxPoints x nSamples array with randomly
%                           sampled POC flux values
%      paramMatrixDepth   - nFluxPoints depths
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 23 Oct 2024 
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

% Number of depths at which POC flux will be extracted
nFluxPoints = sum(~isnan(profileDepths));

paramMatrixDepth = profileDepths(1:nFluxPoints);
paramMatrixFlux = NaN(nFluxPoints,nSamples);
for iFluxPoint = 1:nFluxPoints
    if (~isnan(profileFlux_mean(iFluxPoint)))
        paramMatrixFlux(iFluxPoint,:) = normrnd(profileFlux_mean(iFluxPoint),...
            profileFlux_error(iFluxPoint),[1 nSamples]);
        % Replace negative values with the min value larger than 0 in the
        % string of values generated for that particular flux point
        thisRow = paramMatrixFlux(iFluxPoint,:);
        
        if numel(thisRow) > 1
            % Only replace negative values if there are positive values in the row
            if any(thisRow > 0)
                thisRow(thisRow < 0) = min(thisRow(thisRow > 0));
            end
        else
            % For a single element, only replace it if it is negative and there's a positive value
            if thisRow < 0
                thisRow = 0;
            end
        end

        paramMatrixFlux(iFluxPoint,:) = thisRow;
    end
end

end % samplePocFluxWithLH
