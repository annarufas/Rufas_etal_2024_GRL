function [classic] = calculateBcpMetricsFromTrapAndRadionuclide(...
    isMeansOfMeans,isLogTransformed,isFluxNormalised,choiceZref)

% CALCULATEBCPMETRICSFROMTRAPANDRADIONUCLIDE Computes three key metrics of the 
% biological carbon pump (BCP) —transfer efficiency (Teff), Martin's b, 
% remineralisation depth scale (z*)— using the sediment trap and radionuclide 
% particulate organic carbon (POC) flux measurements compiled for this study. 
% To facilitate smoother curve fitting for the parameters b and z*, the function 
% bins the data into regular 5-meter depth intervals.
%
%    INPUT:
%       isMeansOfMeans   - Boolean indicating whether to calculate the metrics as the 
%                          mean of monthly averages or from annual POC flux values.
%       isLogTransformed - Boolean indicating whether to log-transform depth and POC 
%                          flux data before fitting b and z*.
%       isFluxNormalised - Boolean indicating whether to normalise POC flux
%                          to the value at zref.
%       choiceZref       - Boolean indicating the reference depth choice, 
%                          where 1=closest value to 100, 2=zeu, 3=inflexion point.
%
%    OUTPUT:
%       classic - local annual average values and standard
%                 deviations of b, z* and Teff 100 to 1000 m
%
%   This script uses these external functions: 
%       propagateErrorWithMCforTeff.m            - custom function
%       propagateErrorWithMCforMartinbAndZstar.m - custom function
%       binPocFluxDataAtRegularDepthIntervals.m  - custom function
%       extractDataBelowZref.m                   - custom function
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

%% Presets

% Filename declarations 
filenameInputMonthlyPocFlux = 'pocflux_compilation_monthly.mat';
filenameInputTimeseriesInformation = 'timeseries_station_information.mat';
filenameOutputBcpMetrics = 'bcpmetrics_classic.mat';

% Define parameters
MAX_NUM_DEPTHS_PER_PROFILE = 100;
MOLAR_MASS_CARBON = 12.011; % g mol-1

% Set the seed for the random number generator
rng(0); % this will create a repeatble sequence of random numbers every time randn or randi are called

% Load the trap and radionuclide compilation (mmol C m-2 d-1)
load(fullfile('.','data','processed',filenameInputMonthlyPocFlux),...
    'classicMonthlyDhAvg','classicMonthlyDhErrTot','classicMonthlyProfileAvg',...
    'classicMonthlyProfileErrTot','classicMonthlyProfileDepths',...
    'classicMonthlyProfileN','classicAnnualDhAvg','classicAnnualDhErrTot',...
    'classicAnnualProfileAvg','classicAnnualProfileErrTot','classicAnnualProfileN',...
    'classicAnnualProfileDepths')

% Load information on stations
load(fullfile('.','data','processed',filenameInputTimeseriesInformation),...
    'NUM_LOCS','LOC_LATS','LOC_LONS')

%% Euphotic layer depth (needed to calculate zref)

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

if (~isMeansOfMeans) % calculate annual mean values
    qZeuAnnual = mean(qZeuMonthly,2,'omitnan');
end

%% Calculate Teff

if isMeansOfMeans
    
    % Offload POC flux data to be passed to the function computing Teff 
    arrayFlux = NaN(2,12,NUM_LOCS,2); % mmol C m-2 d-1, 1st dim: 1=base zeu, 2=base zmeso / 4th dim: 1=avg, 2=err
    arrayFlux(1,:,:,1) = squeeze(classicMonthlyDhAvg(1,:,:)); % z x 12 x NUM_LOCS
    arrayFlux(2,:,:,1) = squeeze(classicMonthlyDhAvg(2,:,:)); 
    arrayFlux(1,:,:,2) = squeeze(classicMonthlyDhErrTot(1,:,:));
    arrayFlux(2,:,:,2) = squeeze(classicMonthlyDhErrTot(2,:,:)); 

    fprintf('\nInitiate calculation of Teff...\n')
    [teffAnnual] = propagateErrorWithMCforTeff(arrayFlux,isMeansOfMeans);
    fprintf('\n...done.\n')

    classic.teff.ave      = teffAnnual(:,1);
    classic.teff.stdevupp = teffAnnual(:,2);
    classic.teff.stdevlow = teffAnnual(:,3);
    classic.teff.max      = teffAnnual(:,4);
    classic.teff.min      = teffAnnual(:,5);
    
else
    
    % Offload POC flux data to be passed to the function computing Teff 
    arrayFlux = NaN(2,NUM_LOCS,2); % mmol C m-2 d-1, 1st dim: 1=base zeu, 2=base zmeso / 3rd dim: 1=avg, 2=err
    arrayFlux(1,:,1) = squeeze(classicAnnualDhAvg(1,:)); 
    arrayFlux(2,:,1) = squeeze(classicAnnualDhAvg(2,:));
    arrayFlux(1,:,2) = squeeze(classicAnnualDhErrTot(1,:)); 
    arrayFlux(2,:,2) = squeeze(classicAnnualDhErrTot(2,:));

    fprintf('\nInitiate calculation of Teff...\n')       
    [teffAnnual] = propagateErrorWithMCforTeff(arrayFlux,isMeansOfMeans);
    fprintf('\n...done.\n')

    classic.teff.ave      = teffAnnual(:,1);
    classic.teff.stdevupp = teffAnnual(:,2);
    classic.teff.stdevlow = teffAnnual(:,3);
    classic.teff.max      = teffAnnual(:,4);
    classic.teff.min      = teffAnnual(:,5);
    
end

%% Calculate Martin's b and z*

if isMeansOfMeans

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

    fprintf('\nInitiate calculation of b and z*...\n')
    [martinbAnnual,zstarAnnual,martinb_gof,zstar_gof] =...
        propagateErrorWithMCforMartinbAndZstar(arrayDepths,arrayFlux,...
            isMeansOfMeans,isLogTransformed,isFluxNormalised,choiceZref);
    fprintf('\n...done.\n')

else

    % Bin data at 5 m depth regular intervals to smooth curve fitting
    fluxBinned   = NaN(MAX_NUM_DEPTHS_PER_PROFILE,NUM_LOCS,2); % mmol C m-2 d-1, 4th dim: 1=avg, 2=err
    depthsBinned = NaN(MAX_NUM_DEPTHS_PER_PROFILE,NUM_LOCS);
  
    for iLoc = 1:NUM_LOCS 
        [depthsBinned(:,iLoc),fluxBinned(:,iLoc,:)] =... 
            binPocFluxDataAtRegularDepthIntervals(...
                classicAnnualProfileDepths(:,iLoc),...
                classicAnnualProfileAvg(:,iLoc),...
                classicAnnualProfileErrTot(:,iLoc),...
                classicAnnualProfileN(:,iLoc),...
                MAX_NUM_DEPTHS_PER_PROFILE);
    end 
    
    % Extract data from zref and below, to be passed to the function 
    % calculating b and z*
    arrayFlux   = NaN(MAX_NUM_DEPTHS_PER_PROFILE,NUM_LOCS,2); % mg C m-2 d-1, 4th dim: 1=mean, 2=err 
    arrayDepths = NaN(MAX_NUM_DEPTHS_PER_PROFILE,NUM_LOCS);
    
    for iLoc = 1:NUM_LOCS

        [selectedDepths,selectedFluxes,selectedErrors] = ...
            extractDataBelowZref(depthsBinned(:,iLoc),fluxBinned(:,iLoc,1),...
                fluxBinned(:,iLoc,2),choiceZref,qZeuAnnual(iLoc));
  
        nEntriesSelected = numel(selectedDepths);
        if (~isempty(nEntriesSelected) && nEntriesSelected > 0)
            if isFluxNormalised
                selectedFluxes = selectedFluxes./selectedFluxes(1);
                selectedErrors = selectedErrors./selectedErrors(1);
            end
            arrayFlux(1:nEntriesSelected,iLoc,1) = MOLAR_MASS_CARBON.*selectedFluxes; % mmol C m-2 d-1 --> mg C m-2 d-1
            arrayFlux(1:nEntriesSelected,iLoc,2) = MOLAR_MASS_CARBON.*selectedErrors;
            arrayDepths(1:nEntriesSelected,iLoc) = selectedDepths;
        end

    end 
    
    fprintf('\nInitiate calculation of b and z*...\n')
    [martinbAnnual,zstarAnnual,martinb_gof,zstar_gof] =... 
        propagateErrorWithMCforMartinbAndZstar(arrayDepths,arrayFlux,...
            isMeansOfMeans,isLogTransformed,isFluxNormalised,choiceZref);
    fprintf('\n...done.\n')
 
end

classic.martinb.ave      = martinbAnnual(:,1);
classic.martinb.stdevupp = martinbAnnual(:,2);
classic.martinb.stdevlow = martinbAnnual(:,3);
classic.martinb.max      = martinbAnnual(:,4);
classic.martinb.min      = martinbAnnual(:,5);
classic.martinb.gof      = martinb_gof;

classic.zstar.ave        = zstarAnnual(:,1);
classic.zstar.stdevupp   = zstarAnnual(:,2);
classic.zstar.stdevlow   = zstarAnnual(:,3);
classic.zstar.max        = zstarAnnual(:,4);
classic.zstar.min        = zstarAnnual(:,5);
classic.zstar.gof        = zstar_gof;

%% Save all the metrics

save(fullfile('.','data','processed',filenameOutputBcpMetrics),'classic*')

end
