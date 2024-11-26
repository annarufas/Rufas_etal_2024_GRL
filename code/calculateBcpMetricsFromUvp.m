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
%       extractDataFromZrefToZmeso.m             - custom function
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 23 Nov 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Presets

% Filename declarations 
filenameInputUvpProcessedDataset45sc = 'pocflux_bisson_45sc.mat';
filenameInputTimeseriesInformation   = 'timeseries_station_information.mat';
filenameOutputBcpMetrics             = 'bcpmetrics_uvp.mat';

% Set the seed for the random number generator
rng(0); % this will create a repeatble sequence of random numbers every time randn or randi are called

% Load the UVP5 dataset 
load(fullfile('.','data','processed','UVP5',filenameInputUvpProcessedDataset45sc),...
    'LOC_DEPTH_HORIZONS','ecotaxaDepths','uvpMonthlyFluxProfileAvg','uvpMonthlyFluxProfileErr',...
    'uvpMonthlyFluxDhAvg','uvpMonthlyFluxDhErr','uvpAnnualFluxDhAvg','uvpAnnualFluxDhErr',...
    'uvpAnnualFluxProfileAvg','uvpAnnualFluxProfileErr')
                
% Load information on reference depths used to extract trap and
% radionuclide POC flux data for Teff calculation
load(fullfile('.','data','processed',filenameInputTimeseriesInformation),...
    'MAX_ZEU','STATION_NAMES')
nLocs = length(STATION_NAMES);

%% Calculate Teff

if isMeansOfMeans

    % Array to offload POC flux data to be passed to the function computing Teff 
    arrayFlux = NaN(2,12,nLocs,2); % mg C m-2 d-1, 1st dim: 1=base zeu, 2=base zmeso / 4th dim: 1=avg, 2=err
    arrayFlux(1,:,:,1) = squeeze(uvpMonthlyFluxDhAvg(1,:,:)); 
    arrayFlux(2,:,:,1) = squeeze(uvpMonthlyFluxDhAvg(2,:,:)); 
    arrayFlux(1,:,:,2) = squeeze(uvpMonthlyFluxDhErr(1,:,:));  
    arrayFlux(2,:,:,2) = squeeze(uvpMonthlyFluxDhErr(2,:,:)); 

    fprintf('\nInitiate calculation of Teff...\n')
    [teffAnnual,teffMonthly] = propagateErrorWithMCforTeff(arrayFlux,isMeansOfMeans);
    fprintf('\n...done.\n')

    uvp.teff.monthlyave      = teffMonthly(:,:,1);
    uvp.teff.monthlystdevupp = teffMonthly(:,:,2);
    uvp.teff.monthlystdevlow = teffMonthly(:,:,3);
    
else
        
    % Offload POC flux data to be passed to the function computing Teff 
    arrayFlux = NaN(2,nLocs,2); % mg C m-2 d-1, 1st dim: 1=base zeu, 2=base zmeso / 3rd dim: 1=avg, 2=err
    arrayFlux(1,:,1) = squeeze(uvpAnnualFluxDhAvg(1,:)); 
    arrayFlux(2,:,1) = squeeze(uvpAnnualFluxDhAvg(2,:));
    arrayFlux(1,:,2) = squeeze(uvpAnnualFluxDhErr(1,:)); 
    arrayFlux(2,:,2) = squeeze(uvpAnnualFluxDhErr(2,:));

    fprintf('\nInitiate calculation of Teff...\n')       
    [teffAnnual,~] = propagateErrorWithMCforTeff(arrayFlux,isMeansOfMeans);
    fprintf('\n...done.\n')

end

uvp.teff.ave      = teffAnnual(:,1);
uvp.teff.stdevupp = teffAnnual(:,2);
uvp.teff.stdevlow = teffAnnual(:,3);
uvp.teff.max      = teffAnnual(:,4);
uvp.teff.min      = teffAnnual(:,5);

%% Calculate Martin's b and z*

% To smooth curve fitting, reduce resolution of UVP profiles
curveFitDepths = [22.5:50:2000.5]';
        
if isMeansOfMeans

    % Extract data from zref and below, to be passed to the function 
    % calculating b and z*
    arrayFlux   = NaN(numel(ecotaxaDepths),12,nLocs,2); % mg C m-2 d-1, 4th dim: 1=avg, 2=err 
    arrayDepths = NaN(numel(ecotaxaDepths),12,nLocs);
    
    for iLoc = 1:nLocs
        for iMonth = 1:12

            [selectedDepths,selectedFluxes,selectedErrors] =...
                extractDataFromZrefToZmeso(ecotaxaDepths,...
                    uvpMonthlyFluxProfileAvg(:,iMonth,iLoc),...
                    uvpMonthlyFluxProfileErr(:,iMonth,iLoc),...
                    choiceZref,squeeze(LOC_DEPTH_HORIZONS(iMonth,iLoc,:,:)),MAX_ZEU);
            
            % Reduce UVP profiles to target depths only to smooth curve fitting
            [isTarget,~] = ismember(selectedDepths,curveFitDepths);
            selectedFluxes = selectedFluxes(isTarget);  
            selectedErrors = selectedErrors(isTarget);
            selectedDepths = selectedDepths(isTarget); 
            
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

    fprintf('\nInitiate calculation of b and z*...\n')
    [martinbAnnual,zstarAnnual,martinb_gof,zstar_gof,martinbMonthly,zstarMonthly] =... 
        propagateErrorWithMCforMartinbAndZstar(arrayDepths,arrayFlux,...
            isMeansOfMeans,isLogTransformed,isFluxNormalised,choiceZref);
    fprintf('\n...done.\n')
    
    uvp.martinb.monthlyave      = martinbMonthly(:,:,1);
    uvp.martinb.monthlystdevupp = martinbMonthly(:,:,2);
    uvp.martinb.monthlystdevlow = martinbMonthly(:,:,3);
    
    uvp.zstar.monthlyave        = zstarMonthly(:,:,1);
    uvp.zstar.monthlystdevupp   = zstarMonthly(:,:,2);
    uvp.zstar.monthlystdevlow   = zstarMonthly(:,:,3);
  
else

    % Extract data from zref and below, to be passed to the function 
    % calculating b and z*
    arrayFlux   = NaN(MAX_NUM_DEPTHS_PER_PROFILE,nLocs,2); % mg C m-2 d-1, 4th dim: 1=mean, 2=err 
    arrayDepths = NaN(MAX_NUM_DEPTHS_PER_PROFILE,nLocs);
         
    for iLoc = 1:nLocs
        depthBounds = [max(LOC_DEPTH_HORIZONS(:,iLoc,1,1)), max(LOC_DEPTH_HORIZONS(:,iLoc,1,2)), max(LOC_DEPTH_HORIZONS(:,iLoc,1,3));
                       max(LOC_DEPTH_HORIZONS(:,iLoc,2,1)), max(LOC_DEPTH_HORIZONS(:,iLoc,2,2)), max(LOC_DEPTH_HORIZONS(:,iLoc,2,3))];

        [selectedDepths,selectedFluxes,selectedErrors] = ...
            extractDataFromZrefToZmeso(ecotaxaDepths,uvpAnnualFluxProfileAvg(:,iLoc),...
                uvpAnnualFluxProfileErr(:,iLoc),choiceZref,depthBounds,MAX_ZEU);
            
        % Reduce UVP profiles to target depths only to smooth curve fitting
        [isTarget,~] = ismember(selectedDepths,curveFitDepths);
        selectedFluxes = selectedFluxes(isTarget);  
        selectedErrors = selectedErrors(isTarget);
        selectedDepths = selectedDepths(isTarget); 
  
        nEntriesSelected = numel(selectedDepths);
        if (~isempty(nEntriesSelected) && nEntriesSelected > 0)
            if isFluxNormalised
                selectedFluxes = selectedFluxes./selectedFluxes(1);
                selectedErrors = selectedErrors./selectedErrors(1);
            end
            arrayFlux(1:nEntriesSelected,iLoc,1) = selectedFluxes;
            arrayFlux(1:nEntriesSelected,iLoc,2) = selectedErrors;
            arrayDepths(1:nEntriesSelected,iLoc) = selectedDepths;
        end

    end 
    
    fprintf('\nInitiate calculation of b and z*...\n')
    [martinbAnnual,zstarAnnual,martinb_gof,zstar_gof,~,~] =... 
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