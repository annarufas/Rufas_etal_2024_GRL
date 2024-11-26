function [classic] = calculateBcpMetricsFromTrapAndRadCompilation(...
    isMeansOfMeans,isLogTransformed,isFluxNormalised,choiceZref)

% CALCULATEBCPMETRICSFROMTRAPANDRADCOMPILATION Computes three key metrics 
% of the ocean's biological carbon pump (BCP) —transfer efficiency (Teff), 
% Martin's b coefficient, remineralisation depth scale coefficient (z*)— 
% using the sediment trap and radionuclide particulate organic carbon (POC) 
% flux measurements compiled for this study. To facilitate smoother curve 
% fitting for the parameters b and z*, the function bins the data into 
% regular 5-metre depth intervals.
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
%                 deviations of b, z* and Teff 100 to 1000 m, as well as
%                 monthly values
%
%   This script uses these external functions: 
%       propagateErrorWithMCforTeff.m            - custom function
%       propagateErrorWithMCforMartinbAndZstar.m - custom function
%       binPocFluxDataAtRegularDepthIntervals.m  - custom function
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
filenameInputFluxCompilation       = 'pocflux_compilation.mat';
filenameInputTimeseriesInformation = 'timeseries_station_information.mat';
filenameOutputBcpMetrics           = 'bcpmetrics_classic.mat';

% Set the seed for the random number generator
rng(0); % this will create a repeatble sequence of random numbers every time randn or randi are called

% Load the trap and radionuclide compilation (mg C m-2 d-1)
load(fullfile('.','data','processed',filenameInputFluxCompilation),...
    'classicMonthlyDhAvg','classicMonthlyDhErrTot','classicMonthlyProfileAvg',...
    'classicMonthlyProfileErrTot','classicMonthlyProfileDepths',...
    'classicMonthlyProfileN','classicAnnualDhAvg','classicAnnualDhErrTot',...
    'classicAnnualProfileAvg','classicAnnualProfileErrTot','classicAnnualProfileN',...
    'classicAnnualProfileDepths')

% Load information on stations
load(fullfile('.','data','processed',filenameInputTimeseriesInformation),...
    'qZeuMonthly','MAX_NUM_DEPTHS_PER_PROFILE','LOC_DEPTH_HORIZONS','MAX_ZEU')
nLocs = size(qZeuMonthly,2);

%% Calculate Teff

if isMeansOfMeans
    
    % Offload POC flux data to be passed to the function computing Teff 
    arrayFlux = NaN(2,12,nLocs,2); % 1st dim: 1=base zeu, 2=base zmeso / 4th dim: 1=avg, 2=err
    arrayFlux(1,:,:,1) = squeeze(classicMonthlyDhAvg(1,:,:)); % z x 12 x nLocs
    arrayFlux(2,:,:,1) = squeeze(classicMonthlyDhAvg(2,:,:)); 
    arrayFlux(1,:,:,2) = squeeze(classicMonthlyDhErrTot(1,:,:));
    arrayFlux(2,:,:,2) = squeeze(classicMonthlyDhErrTot(2,:,:)); 

    fprintf('\nInitiate calculation of Teff...\n')
    [teffAnnual,teffMonthly] = propagateErrorWithMCforTeff(arrayFlux,isMeansOfMeans);
    fprintf('\n...done.\n')
    
    classic.teff.monthlyave      = teffMonthly(:,:,1);
    classic.teff.monthlystdevupp = teffMonthly(:,:,2);
    classic.teff.monthlystdevlow = teffMonthly(:,:,3);

else
    
    % Offload POC flux data to be passed to the function computing Teff 
    arrayFlux = NaN(2,nLocs,2); % 1st dim: 1=base zeu, 2=base zmeso / 3rd dim: 1=avg, 2=err
    arrayFlux(1,:,1) = squeeze(classicAnnualDhAvg(1,:)); 
    arrayFlux(2,:,1) = squeeze(classicAnnualDhAvg(2,:));
    arrayFlux(1,:,2) = squeeze(classicAnnualDhErrTot(1,:)); 
    arrayFlux(2,:,2) = squeeze(classicAnnualDhErrTot(2,:));

    fprintf('\nInitiate calculation of Teff...\n')       
    [teffAnnual,~] = propagateErrorWithMCforTeff(arrayFlux,isMeansOfMeans);
    fprintf('\n...done.\n')

end

classic.teff.ave      = teffAnnual(:,1);
classic.teff.stdevupp = teffAnnual(:,2);
classic.teff.stdevlow = teffAnnual(:,3);
classic.teff.max      = teffAnnual(:,4);
classic.teff.min      = teffAnnual(:,5);

%% Calculate Martin's b and z*

if isMeansOfMeans

    % Bin data at 5 m depth regular intervals to smooth curve fitting
    fluxBinned   = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,nLocs,2); % 4th dim: 1=avg, 2=err
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
  
    % Extract data from zref to zmeso, to be passed to the function 
    % calculating b and z*
    arrayFlux   = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,nLocs,2); % 4th dim: 1=avg, 2=err 
    arrayDepths = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,nLocs);
    
    for iLoc = 1:nLocs
        for iMonth = 1:12

            [selectedDepths,selectedFluxes,selectedErrors] = ...
                extractDataFromZrefToZmeso(depthsBinned(:,iMonth,iLoc),...
                    fluxBinned(:,iMonth,iLoc,1),fluxBinned(:,iMonth,iLoc,2),...
                    choiceZref,squeeze(LOC_DEPTH_HORIZONS(iMonth,iLoc,:,:)),MAX_ZEU);
             
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
    
    classic.martinb.monthlyave      = martinbMonthly(:,:,1);
    classic.martinb.monthlystdevupp = martinbMonthly(:,:,2);
    classic.martinb.monthlystdevlow = martinbMonthly(:,:,3);
    
    classic.zstar.monthlyave        = zstarMonthly(:,:,1);
    classic.zstar.monthlystdevupp   = zstarMonthly(:,:,2);
    classic.zstar.monthlystdevlow   = zstarMonthly(:,:,3);

else
   
    % Bin data at 5 m depth regular intervals to smooth curve fitting
    fluxBinned   = NaN(MAX_NUM_DEPTHS_PER_PROFILE,nLocs,2); % 4th dim: 1=avg, 2=err
    depthsBinned = NaN(MAX_NUM_DEPTHS_PER_PROFILE,nLocs);
  
    for iLoc = 1:nLocs 
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
    arrayFlux   = NaN(MAX_NUM_DEPTHS_PER_PROFILE,nLocs,2); % 4th dim: 1=mean, 2=err 
    arrayDepths = NaN(MAX_NUM_DEPTHS_PER_PROFILE,nLocs);
    
    for iLoc = 1:nLocs

        depthBounds = [max(LOC_DEPTH_HORIZONS(:,iLoc,1,1)), max(LOC_DEPTH_HORIZONS(:,iLoc,1,2)), max(LOC_DEPTH_HORIZONS(:,iLoc,1,3));
                       max(LOC_DEPTH_HORIZONS(:,iLoc,2,1)), max(LOC_DEPTH_HORIZONS(:,iLoc,2,2)), max(LOC_DEPTH_HORIZONS(:,iLoc,2,3))];

        [selectedDepths,selectedFluxes,selectedErrors] = ...
            extractDataFromZrefToZmeso(depthsBinned(:,iLoc),fluxBinned(:,iLoc,1),...
                fluxBinned(:,iLoc,2),choiceZref,depthBounds,MAX_ZEU);

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
