
% ======================================================================= %
%                                                                         %
% This script investigates the effct of sample size in the calculation of % 
% b and z* from POC flux data from our sediment trap/radionuclide         %
% compialtion. The script has 4 sections:                                 %
%   Section 1 - Presets.                                                  %
%   Section 2 - Load euphotic layer depth (needed to calculate zref).     %
%   Section 3 - Rearrange POC flux data for calculations.                 %
%   Section 4 - Calculate the coefficients b and z* and their variation.  %
%   Section 5 - Save the calculated data.                                 %
%   Section 6 - Calculate relative error.                                 %
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
%   Version 1.0 - Completed 7 Nov 2024                                    %
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
filenameInputMonthlyPocFlux = 'pocflux_compilation_monthly.mat';
filenameInputTimeseriesInformation = 'timeseries_station_information.mat';
filenameOutputSampleSizeStats = 'sampleSizeStats.mat';

% Load station information
load(fullfile('.','data','processed',filenameInputTimeseriesInformation),...
    'NUM_LOCS','LOC_LATS','LOC_LONS')

% Load the trap and radionuclide compilation (mmol C m-2 d-1)
load(fullfile('.','data','processed',filenameInputMonthlyPocFlux),...
    'classicMonthlyDhAvg','classicMonthlyDhErrTot','classicMonthlyDhN',...
    'classicMonthlyProfileAvg','classicMonthlyProfileErrTot','classicMonthlyProfileN',...    
    'classicMonthlyProfileDepths')

% Define parameters for sampling
NUM_SIMULATED_PROFILES = 50; 
NUM_SUBSAMPLING = 100;

% Define other parameters
MAX_NUM_DEPTHS_PER_PROFILE = 100;
MOLAR_MASS_CARBON = 12.011; % g mol-1

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - LOAD EUPHOTIC LAYER DEPTH (NEEDED TO CALCULATE ZREF)
% -------------------------------------------------------------------------

% Load the global-ocean euphotic layer depth product calculated from 
% CMEMS kd and extract data at our locations
load(fullfile('.','data','interim','zeu_calculated_onepercentpar0.mat'),'zeu','zeu_lat','zeu_lon')

% Query points for interpolation
qLats = LOC_LATS;
qLons = LOC_LONS;

% Original data grid
[X,Y,T] = ndgrid(zeu_lat,zeu_lon,(1:12)');

% Interpolant 
F = griddedInterpolant(X, Y, T, zeu, 'linear'); 

% Extract data for our locations defined by qLats and qLons
qZeuMonthly = NaN(length(qLats),12);
for iLoc = 1:length(qLats)
    [qX,qY,qT] = ndgrid(qLats(iLoc),qLons(iLoc),(1:12)');
    qZeuMonthly(iLoc,:) = F(qX,qY,qT);
end
qZeuAnnual = mean(qZeuMonthly,2,'omitnan');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - REARRANGE POC FLUX DATA FOR CALCULATIONS
% -------------------------------------------------------------------------

% Bin data at 5 m depth regular intervals to smooth curve fitting
fluxBinned   = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,NUM_LOCS,2); % mmol C m-2 d-1, 4th dim: 1=avg, 2=err
depthsBinned = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,NUM_LOCS);

for iLoc = 1:NUM_LOCS
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
arrayFlux   = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,NUM_LOCS,2); % mg C m-2 d-1, 4th dim: 1=avg, 2=err 
arrayDepths = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,NUM_LOCS);

for iLoc = 1:NUM_LOCS
    for iMonth = 1:12

        [selectedDepths,selectedFluxes,selectedErrors] = ...
            extractDataBelowZref(...
                depthsBinned(:,iMonth,iLoc),fluxBinned(:,iMonth,iLoc,1),...
                fluxBinned(:,iMonth,iLoc,2),choiceZref,qZeuMonthly(iLoc,iMonth));

        nEntriesSelected = numel(selectedDepths);
        if (~isempty(nEntriesSelected) && nEntriesSelected > 0) 
            if isFluxNormalised
                selectedFluxes = selectedFluxes./selectedFluxes(1);
                selectedErrors = selectedErrors./selectedErrors(1);
            end
            arrayFlux(1:nEntriesSelected,iMonth,iLoc,1) = MOLAR_MASS_CARBON.*selectedFluxes; % mmol C m-2 d-1 --> mg C m-2 d-1
            arrayFlux(1:nEntriesSelected,iMonth,iLoc,2) = MOLAR_MASS_CARBON.*selectedErrors;
            arrayDepths(1:nEntriesSelected,iMonth,iLoc) = selectedDepths;
        end

    end
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - CALCULATE THE COEFFICIENTS B AND Z* AND THEIR VARIATION 
% -------------------------------------------------------------------------

% Initialise output arrays
martinbSimulatedProfile = NaN(NUM_LOCS,NUM_SUBSAMPLING,NUM_SIMULATED_PROFILES,2); % 1=mean, 2=GOF
zstarSimulatedProfile   = NaN(NUM_LOCS,NUM_SUBSAMPLING,NUM_SIMULATED_PROFILES,2); % 1=mean, 2=GOF
martinbOverall          = NaN(NUM_LOCS,2); % 1=mean, 2=std
zstarOverall            = NaN(NUM_LOCS,2); % 1=mean, 2=std

for iLoc = 1:NUM_LOCS
    fprintf('\nLocation %d',iLoc)
    
    for iSubsample = 1:NUM_SUBSAMPLING
        
        % Randomly choose 50 times a number between 1 and 12 (since we've 
        % got 12 "mean" profiles (=12 months))
        randomSamples = randi([1, 12], 1, NUM_SIMULATED_PROFILES)';

        for iSimulatedProfile = 1:NUM_SIMULATED_PROFILES
            profileFlux_mean  = arrayFlux(:,randomSamples(iSimulatedProfile),iLoc,1); % mg C m-2 d-1
            profileFlux_error = arrayFlux(:,randomSamples(iSimulatedProfile),iLoc,2); % mg C m-2 d-1
            profileDepths     = arrayDepths(:,randomSamples(iSimulatedProfile),iLoc); % m
            
            % Only proceed if there are flux data and the first depth is <= 200 m
            if (sum(~isnan(profileFlux_mean)) > 0 && profileDepths(1) <= 200)
    
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

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - SAVE THE CALCULATED DATA
% -------------------------------------------------------------------------

save(fullfile('.','data','processed',filenameOutputSampleSizeStats),...
    'martinbOverall','zstarOverall')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 6 - CALCULATE RELATIVE ERROR
% -------------------------------------------------------------------------

martinbRe = 100.*(martinbOverall(:,2)./martinbOverall(:,1)); 
zstarRe   = 100.*(zstarOverall(:,2)./zstarOverall(:,1)); 

fprintf('\nThe relative error in b for our 6 locations is %2.0f, %2.0f, %2.0f, %2.0f, %2.0f, %2.0f percent',...
    martinbRe(1),martinbRe(2),martinbRe(3),martinbRe(4),martinbRe(5),martinbRe(6))
fprintf('\nThe relative error in z* for our 6 locations is %2.0f, %2.0f, %2.0f, %2.0f, %2.0f, %2.0f percent',...
    zstarRe(1),zstarRe(2),zstarRe(3),zstarRe(4),zstarRe(5),zstarRe(6))   
