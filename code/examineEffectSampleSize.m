
% ======================================================================= %
%                                                                         %
% This script investigates the effct of sample size in the calculation of % 
% b and z* from POC flux data from our sediment trap/radionuclide         %
% compialtion. The script has 4 sections:                                 %
%   Section 1 - Presets.                                                  %
%   Section 2 - Rearrange POC flux data for calculations.                 %
%   Section 3 - Calculate the coefficients b and z* and their variation.  %
%   Section 4 - Calculate relative error.                                 %
%                                                                         %
%   This script uses these external functions:                            %
%       binPocFluxDataAtRegularDepthIntervals.m - custom function         %
%       extractDataBelowZref.m                  - custom function         %
%       samplePocFluxWithLH.m                   - custom function         %
%       solveMartinbAndZstar.m                  - custom function         %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 13 Nov 2024                                   %
%                                                                         %
% ======================================================================= %

close all; clear all; clc
addpath(genpath('./data/raw/'));
addpath(genpath('./data/processed/'));
addpath(genpath('./code/'));
addpath(genpath('./resources/internal/'));

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

% Choices to be applied to sediment trap/radionuclide data and UVP5 data
% when fitting b and z* (best choices determined after running
% 'processPocFluxFits.m')
isMeansOfMeans = 1;   % 0=calculate metrics from annual POC flux values, 1=calculate metrics as means of monthly metric means 
isLogTransformed = 0; % 0=do not log-log transform depth and POC flux data before fitting, 1=do transform 
isFluxNormalised = 0; % 0=do not normalise POC flux values to value at zref, 1=do normalise
choiceZref = 2;       % reference depth, 1=closest value to 100, 2=zeu, 3=inflexion point

% Path and filename declarations
filenameInputPocFluxCompilation    = 'pocflux_compilation.mat';
filenameInputTimeseriesInformation = 'timeseries_station_information.mat';
filenameOutputSampleSizeStats      = 'sampleSizeStats.mat';

% Load station information
load(fullfile('.','data','processed',filenameInputTimeseriesInformation),...
    'qZeuMonthly','MAX_NUM_DEPTHS_PER_PROFILE','MAX_ZEU')
nLocs = size(qZeuMonthly,2);

% Load the trap and radionuclide compilation (mmol C m-2 d-1)
load(fullfile('.','data','processed',filenameInputPocFluxCompilation),...
    'classicMonthlyDhAvg','classicMonthlyDhErrTot','classicMonthlyDhN',...
    'classicMonthlyProfileAvg','classicMonthlyProfileErrTot','classicMonthlyProfileN',...    
    'classicMonthlyProfileDepths')

% Define parameters for sampling
NUM_SIMULATED_PROFILES = 50; 
NUM_SUBSAMPLING = 100;

% Set the seed for the random number generator
rng(0); % this will create a repeatble sequence of random numbers every time randn or randi are called

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - REARRANGE POC FLUX DATA FOR CALCULATIONS
% -------------------------------------------------------------------------

% Bin data at 5 m depth regular intervals to smooth curve fitting
fluxBinned   = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,nLocs,2); % mg C m-2 d-1, 4th dim: 1=avg, 2=err
depthsBinned = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,nLocs);

for iLoc = 1:nLocs
    for iMonth = 1:12
        [depthsBinned(:,iMonth,iLoc),fluxBinned(:,iMonth,iLoc,:)] =... 
            binPocFluxDataAtRegularDepthIntervals(...
                classicMonthlyProfileDepths(:,iMonth,iLoc),...
                classicMonthlyProfileAvg(:,iMonth,iLoc),...
                classicMonthlyProfileErrTot(:,iMonth,iLoc),...
                classicMonthlyProfileN(:,iMonth,iLoc),...
                MAX_NUM_DEPTHS_PER_PROFILE);
    end
end

% Extract data from zref and below, to be passed to the function 
% calculating b and z*
arrayFlux   = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,nLocs,2); % mg C m-2 d-1, 4th dim: 1=avg, 2=err 
arrayDepths = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,nLocs);

for iLoc = 1:nLocs
    for iMonth = 1:12

        [selectedDepths,selectedFluxes,selectedErrors] = ...
            extractDataBelowZref(depthsBinned(:,iMonth,iLoc),...
                fluxBinned(:,iMonth,iLoc,1),fluxBinned(:,iMonth,iLoc,2),...
                choiceZref,qZeuMonthly(iMonth,iLoc),MAX_ZEU);

        nEntriesSelected = numel(selectedDepths);
        if (~isempty(nEntriesSelected) && nEntriesSelected > 0) 
            if isFluxNormalised
                selectedFluxes = selectedFluxes./selectedFluxes(1);
                selectedErrors = selectedErrors./selectedErrors(1);
            end
            arrayFlux(1:nEntriesSelected,iMonth,iLoc,1) = selectedFluxes;
            arrayFlux(1:nEntriesSelected,iMonth,iLoc,2) = selectedErrors;
            arrayDepths(1:nEntriesSelected,iMonth,iLoc) = selectedDepths;
        end

    end
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - CALCULATE THE COEFFICIENTS B AND Z* AND THEIR VARIATION 
% -------------------------------------------------------------------------

% Initialise output arrays
martinbSimulatedProfile = NaN(nLocs,NUM_SUBSAMPLING,NUM_SIMULATED_PROFILES,2); % 1=mean, 2=GOF
zstarSimulatedProfile   = NaN(nLocs,NUM_SUBSAMPLING,NUM_SIMULATED_PROFILES,2); % 1=mean, 2=GOF
martinbOverall          = NaN(nLocs,2); % 1=mean, 2=std
zstarOverall            = NaN(nLocs,2); % 1=mean, 2=std

for iLoc = 1:nLocs
    fprintf('\nLocation %d',iLoc)
    
    for iSubsample = 1:NUM_SUBSAMPLING
        
        % Randomly choose 50 times a number between 1 and 12 (since we've 
        % got 12 "mean" profiles (=12 months))
        randomSamples = randi([1, 12], 1, NUM_SIMULATED_PROFILES)';

        for iSimulatedProfile = 1:NUM_SIMULATED_PROFILES
            profileFlux_mean  = arrayFlux(:,randomSamples(iSimulatedProfile),iLoc,1); % mg C m-2 d-1
            profileFlux_error = arrayFlux(:,randomSamples(iSimulatedProfile),iLoc,2); % mg C m-2 d-1
            profileDepths     = arrayDepths(:,randomSamples(iSimulatedProfile),iLoc); % m
            
            % Only proceed if there are flux data
            if (sum(~isnan(profileFlux_mean)))
    
                % Randomly sample POC flux values from a normal distribution
                [f,z] = samplePocFluxWithLH(profileFlux_mean,profileFlux_error,profileDepths,1);

                % Fit b and z*
                [martinModel,zstarModel] = solveMartinbAndZstar(isLogTransformed,f,z);
                martinbSimulatedProfile(iLoc,iSubsample,iSimulatedProfile,1:2) = martinModel(1:2);
                zstarSimulatedProfile(iLoc,iSubsample,iSimulatedProfile,1:2) = zstarModel(1:2);
                
            end

        end % iSimulatedProfile
    end % iSubsample
    
    martinbSimulatedProfile_flattened = squeeze(martinbSimulatedProfile(iLoc,:,:,1));
    zstarSimulatedProfile_flattened = squeeze(zstarSimulatedProfile(iLoc,:,:,1));
    
    martinbOverall(iLoc,1) = mean(martinbSimulatedProfile_flattened(:),'omitnan');
    martinbOverall(iLoc,2) = std(martinbSimulatedProfile_flattened(:),'omitnan');
    zstarOverall(iLoc,1) = mean(zstarSimulatedProfile_flattened(:),'omitnan');
    zstarOverall(iLoc,2) = std(zstarSimulatedProfile_flattened(:),'omitnan');
    
end % iLoc

save(fullfile('.','data','processed',filenameOutputSampleSizeStats),...
    'martinbOverall','zstarOverall')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - CALCULATE RELATIVE ERROR
% -------------------------------------------------------------------------

martinbRe = 100.*(martinbOverall(:,2)./martinbOverall(:,1)); 
zstarRe   = 100.*(zstarOverall(:,2)./zstarOverall(:,1)); 

fprintf('\nThe relative error in b for our 6 locations is %2.0f, %2.0f, %2.0f, %2.0f, %2.0f, %2.0f percent',...
    martinbRe(1),martinbRe(2),martinbRe(3),martinbRe(4),martinbRe(5),martinbRe(6))
fprintf('\nThe relative error in z* for our 6 locations is %2.0f, %2.0f, %2.0f, %2.0f, %2.0f, %2.0f percent',...
    zstarRe(1),zstarRe(2),zstarRe(3),zstarRe(4),zstarRe(5),zstarRe(6))   
