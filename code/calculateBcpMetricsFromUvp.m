function [uvp] = calculateBcpMetricsFromUvp(isMeansOfMeans,isLogTransformed,...
    isFluxNormalised,choiceZref)

% CALCULATEBCPMETRICSFROMUVP Calculates three metrics of the BCP (Teff, b 
% and z*) using the UVP-5 derived POC flux data calculated using the model 
% of Bisson et al. (2022).
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
%   OUTPUT:
%       uvp - local annual average values and standard
%             deviations of b, z* and Teff 100 to 1000 m
%
%   This script uses these external functions: 
%       propagateErrorWithMCforTeff.m            - custom function
%       propagateErrorWithMCforMartinbAndZstar.m - custom function
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 16 Oct 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Presets

% Filename declarations 
filenameInputUvpProcessedDataset45sc = 'pocflux_bisson_45sc_monthly_and_annual_all_depths.mat';
filenameInputTimeseriesInformation = 'timeseries_station_information.mat';
filenameOutputBcpMetrics = 'bcpmetrics_uvp.mat';

% Define parameters
MOLAR_MASS_CARBON = 12.011; % g mol-1

% Set the seed for the random number generator
rng(0); % this will create a repeatble sequence of random numbers every time randn or randi are called

% Load the UVP5 dataset 
load(fullfile('.','data','processed','UVP5',filenameInputUvpProcessedDataset45sc),...
    'ecotaxaDepths','uvpMonthlyFluxProfileAvg','uvpMonthlyFluxProfileErr',...
    'uvpMonthlyFluxDhAvg','uvpMonthlyFluxDhErr','uvpAnnualFluxDhAvg','uvpAnnualFluxDhErr',...
    'uvpAnnualFluxProfileAvg','uvpAnnualFluxProfileErr')
                
% Load information on reference depths used to extract trap and
% radionuclide POC flux data for Teff calculation
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

    % Array to offload POC flux data to be passed to the function computing Teff 
    arrayFlux = NaN(2,12,NUM_LOCS,2); % mmol C m-2 d-1, 1st dim: 1=base zeu, 2=base zmeso / 4th dim: 1=avg, 2=err
    arrayFlux(1,:,:,1) = squeeze(uvpMonthlyFluxDhAvg(1,:,:)); 
    arrayFlux(2,:,:,1) = squeeze(uvpMonthlyFluxDhAvg(2,:,:)); 
    arrayFlux(1,:,:,2) = squeeze(uvpMonthlyFluxDhErr(1,:,:));  
    arrayFlux(2,:,:,2) = squeeze(uvpMonthlyFluxDhErr(2,:,:)); 

    fprintf('\nInitiate calculation of Teff...\n')
    [teffAnnual] = propagateErrorWithMCforTeff(arrayFlux,isMeansOfMeans);
    fprintf('\n...done.\n')

    uvp.teff.ave      = teffAnnual(:,1);
    uvp.teff.stdevupp = teffAnnual(:,2);
    uvp.teff.stdevlow = teffAnnual(:,3);
    uvp.teff.max      = teffAnnual(:,4);
    uvp.teff.min      = teffAnnual(:,5);

else
        
    % Offload POC flux data to be passed to the function computing Teff 
    arrayFlux = NaN(2,NUM_LOCS,2); % mmol C m-2 d-1, 1st dim: 1=base zeu, 2=base zmeso / 3rd dim: 1=avg, 2=err
    arrayFlux(1,:,1) = squeeze(uvpAnnualFluxDhAvg(1,:)); 
    arrayFlux(2,:,1) = squeeze(uvpAnnualFluxDhAvg(2,:));
    arrayFlux(1,:,2) = squeeze(uvpAnnualFluxDhErr(1,:)); 
    arrayFlux(2,:,2) = squeeze(uvpAnnualFluxDhErr(2,:));

    fprintf('\nInitiate calculation of Teff...\n')       
    [teffAnnual] = propagateErrorWithMCforTeff(arrayFlux,isMeansOfMeans);
    fprintf('\n...done.\n')

    uvp.teff.ave      = teffAnnual(:,1);
    uvp.teff.stdevupp = teffAnnual(:,2);
    uvp.teff.stdevlow = teffAnnual(:,3);
    uvp.teff.max      = teffAnnual(:,4);
    uvp.teff.min      = teffAnnual(:,5);
    
end


%% Calculate Martin's b and z*

% Reduce UVP profiles to target depths only to smooth curve fitting
metricCurveFitDepths = [25,50,75,100,150,200,250,300,500,750,1000,1500,2000]'; % to calculate b and z*
izFitDepths = NaN(numel(metricCurveFitDepths),1);
for iDepth = 1:numel(metricCurveFitDepths)
    [~,izFitDepths(iDepth)] = min(abs(ecotaxaDepths - metricCurveFitDepths(iDepth)));
end 

nCurveFitDepths = numel(metricCurveFitDepths); 

if isMeansOfMeans
    
    % Flux for curve fitting
    fluxCurveFit = NaN(nCurveFitDepths,12,NUM_LOCS,2); % mmol C m-2 d-1, 4th dim: 1=avg, 2=err
    fluxCurveFit(:,:,:,1) = squeeze(uvpMonthlyFluxProfileAvg(izFitDepths,:,:));  
    fluxCurveFit(:,:,:,2) = squeeze(uvpMonthlyFluxProfileErr(izFitDepths,:,:));
    depthsCurveFit = ecotaxaDepths(izFitDepths); 

    % Extract data from zref and below, to be passed to the function 
    % calculating b and z*
    arrayFlux   = NaN(nCurveFitDepths,12,NUM_LOCS,2); % mg C m-2 d-1, 4th dim: 1=avg, 2=err 
    arrayDepths = NaN(nCurveFitDepths,12,NUM_LOCS);
    
    for iLoc = 1:NUM_LOCS
        for iMonth = 1:12

            [selectedDepths,selectedFluxes,selectedErrors] =...
                extractDataBelowZref(depthsCurveFit,fluxCurveFit(:,iMonth,iLoc,1),...
                    fluxCurveFit(:,iMonth,iLoc,2),choiceZref,qZeuMonthly(iLoc,iMonth));
             
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
    
    % Flux and depths for curve fitting
    fluxCurveFit = NaN(nCurveFitDepths,NUM_LOCS,2); % mmol C m-2 d-1, 4th dim: 1=avg, 2=err
    fluxCurveFit(:,:,1) = squeeze(uvpAnnualFluxProfileAvg(izFitDepths,:));  
    fluxCurveFit(:,:,2) = squeeze(uvpAnnualFluxProfileErr(izFitDepths,:));
    depthsCurveFit = ecotaxaDepths(izFitDepths); 
    
    % Extract data from zref and below, to be passed to the function 
    % calculating b and z*
    arrayFlux   = NaN(MAX_NUM_DEPTHS_PER_PROFILE,NUM_LOCS,2); % mg C m-2 d-1, 4th dim: 1=mean, 2=err 
    arrayDepths = NaN(MAX_NUM_DEPTHS_PER_PROFILE,NUM_LOCS);
    
    for iLoc = 1:NUM_LOCS

        [selectedDepths,selectedFluxes,selectedErrors] = ...
            extractDataBelowZref(depthsCurveFit,fluxCurveFit(:,iLoc,1),...
                fluxCurveFit(:,iLoc,2),choiceZref,qZeuAnnual(iLoc));
  
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

uvp.martinb.ave      = martinbAnnual(:,1);
uvp.martinb.stdevupp = martinbAnnual(:,2);
uvp.martinb.stdevlow = martinbAnnual(:,3);
uvp.martinb.max      = martinbAnnual(:,4);
uvp.martinb.min      = martinbAnnual(:,5);
uvp.martinb.gof      = martinb_gof;

uvp.zstar.ave        = zstarAnnual(:,1);
uvp.zstar.stdevupp   = zstarAnnual(:,2);
uvp.zstar.stdevlow   = zstarAnnual(:,3);
uvp.zstar.max        = zstarAnnual(:,4);
uvp.zstar.min        = zstarAnnual(:,5);
uvp.zstar.gof        = zstar_gof;

%% Save metrics

save(fullfile('.','data','processed',filenameOutputBcpMetrics),'uvp*')

end