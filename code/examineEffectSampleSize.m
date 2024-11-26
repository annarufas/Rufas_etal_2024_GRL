
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
%       extractDataFromZrefToZmeso.m            - custom function         %
%       samplePocFluxWithLH.m                   - custom function         %
%       solveMartinbAndZstar.m                  - custom function         %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 23 Nov 2024                                   %
%                                                                         %
% ======================================================================= %

close all; clear all; clc
addpath(genpath('./data/raw/'));
addpath(genpath('./data/processed/'));
addpath(genpath('./code/'));
addpath(genpath('./resources/internal/'));
addpath(genpath('./resources/external/'));

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
filenameInputMetricsData           = 'bcpmetrics_all.mat';

% Load station information
load(fullfile('.','data','processed',filenameInputTimeseriesInformation),...
    'qZeuMonthly','MAX_NUM_DEPTHS_PER_PROFILE','LOC_DEPTH_HORIZONS','MAX_ZEU','STATION_NAMES')
nLocs = length(STATION_NAMES);

% Load the trap and radionuclide compilation (mmol C m-2 d-1)
load(fullfile('.','data','processed',filenameInputPocFluxCompilation),...
    'classicMonthlyDhAvg','classicMonthlyDhErrTot','classicMonthlyDhN',...
    'classicMonthlyProfileAvg','classicMonthlyProfileErrTot','classicMonthlyProfileN',...    
    'classicMonthlyProfileDepths')

% Load the metrics array
load(fullfile('.','data','processed',filenameInputMetricsData),'martinbMonthlyTraAndRad')

% Define parameters for sampling
NUM_SIMULATED_PROFILES = 50; 
NUM_SUBSAMPLING = 1000;

% Set the seed for the random number generator
rng(0); % this will create a repeatble sequence of random numbers every time randn or randi are called

% Indexes to locations
iE = 1; % EqPac
iO = 2; % OSP
iP = 3; % PAP-SO
iB = 4; % BATS/OFP
iHo = 5; % HOT/ALOHA
iHa = 6; % HAUSGARTEN

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

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - CALCULATE THE COEFFICIENTS B AND Z* AND THEIR VARIATION 
% -------------------------------------------------------------------------

% Initialise output arrays
martinbSimulatedProfile = NaN(nLocs,NUM_SIMULATED_PROFILES,NUM_SUBSAMPLING,2); % 1=mean, 2=GOF
zstarSimulatedProfile   = NaN(nLocs,NUM_SIMULATED_PROFILES,NUM_SUBSAMPLING,2); % 1=mean, 2=GOF
martinbSimulatedOverall = NaN(nLocs,2); % 1=mean, 2=std
zstarSimulatedOverall   = NaN(nLocs,2); % 1=mean, 2=std
idxMonth                = NaN(nLocs,NUM_SIMULATED_PROFILES);

for iLoc = 1:nLocs
    fprintf('\nLocation %d',iLoc)

    % Randomly choose 50 times a number between 1 and 12 (since we've 
    % got 12 "mean" profiles (=12 months))
    localData = squeeze(arrayFlux(:,:,iLoc,1));
    nanColumns = all(isnan(localData), 1);
    idxValidColumns = find(~nanColumns);
    
     % Since there are fewer than 50 non-NaN columns, sample with replacement
    randomIndices = idxValidColumns(randi(numel(idxValidColumns), 1, NUM_SIMULATED_PROFILES, true));
    idxMonth(iLoc,:) = randomIndices;
    
    for iSimulatedProfile = 1:NUM_SIMULATED_PROFILES
        profileFlux_mean  = arrayFlux(:,randomIndices(iSimulatedProfile),iLoc,1); % mg C m-2 d-1
        profileFlux_error = arrayFlux(:,randomIndices(iSimulatedProfile),iLoc,2); % mg C m-2 d-1
        profileDepths     = arrayDepths(:,randomIndices(iSimulatedProfile),iLoc); % m
            
        for iSubsample = 1:NUM_SUBSAMPLING
            if sum(~isnan(profileFlux_mean)) > 0
    
                % Randomly sample POC flux values from a normal distribution
                [f,z] = samplePocFluxWithLH(profileFlux_mean,profileFlux_error,profileDepths,1);

                % Fit b and z*
                [martinModel,zstarModel] = solveMartinbAndZstar(isLogTransformed,f,z);
                martinbSimulatedProfile(iLoc,iSimulatedProfile,iSubsample,1:2) = martinModel(1:2);
                zstarSimulatedProfile(iLoc,iSimulatedProfile,iSubsample,1:2) = zstarModel(1:2);
                
            end

        end % iSubsample 
    end % iSimulatedProfile
    
    martinbSimulatedProfile_flattened = squeeze(martinbSimulatedProfile(iLoc,:,:,1));
    zstarSimulatedProfile_flattened = squeeze(zstarSimulatedProfile(iLoc,:,:,1));
    
    martinbSimulatedOverall(iLoc,1) = mean(martinbSimulatedProfile_flattened(:),'omitnan');
    martinbSimulatedOverall(iLoc,2) = std(martinbSimulatedProfile_flattened(:),'omitnan');
    zstarSimulatedOverall(iLoc,1) = mean(zstarSimulatedProfile_flattened(:),'omitnan');
    zstarSimulatedOverall(iLoc,2) = std(zstarSimulatedProfile_flattened(:),'omitnan');
    
end % iLoc

save(fullfile('.','data','processed',filenameOutputSampleSizeStats),...
    'martinbSimulatedOverall','zstarSimulatedOverall',...
    'martinbSimulatedProfile','zstarSimulatedProfile',...
    'idxMonth')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - CALCULATE RELATIVE ERROR IN SIMULATED B AND CORRELATION
% COEFFICIENTS
% -------------------------------------------------------------------------

load(fullfile('.','data','processed',filenameOutputSampleSizeStats),...
    'martinbSimulatedOverall','zstarSimulatedOverall',...
    'martinbSimulatedProfile','zstarSimulatedProfile',...
    'idxMonth')

% Relative error
martinbRelativeError = 100.*(martinbSimulatedOverall(:,2)./martinbSimulatedOverall(:,1)); 
zstarRelativeError   = 100.*(zstarSimulatedOverall(:,2)./zstarSimulatedOverall(:,1)); 
fprintf('\nThe relative error in b for our 6 locations is %2.0f, %2.0f, %2.0f, %2.0f, %2.0f, %2.0f percent',...
    martinbRelativeError(1),martinbRelativeError(2),martinbRelativeError(3),martinbRelativeError(4),martinbRelativeError(5),martinbRelativeError(6))
fprintf('\nThe relative error in z* for our 6 locations is %2.0f, %2.0f, %2.0f, %2.0f, %2.0f, %2.0f percent',...
    zstarRelativeError(1),zstarRelativeError(2),zstarRelativeError(3),zstarRelativeError(4),zstarRelativeError(5),zstarRelativeError(6))   

% Rescue true value
martinbTrue = NaN(nLocs,NUM_SIMULATED_PROFILES,NUM_SUBSAMPLING);
for iLoc = 1:nLocs
    idxTrueMartinb = repmat(squeeze(idxMonth(iLoc,:))',[1,NUM_SUBSAMPLING]);
    for iRow = 1:size(idxTrueMartinb,1)
        for iCol = 1:size(idxTrueMartinb,2)
            thisMartinb = idxTrueMartinb(iRow,iCol);
            martinbTrue(iLoc,iRow,iCol) = martinbMonthlyTraAndRad(iLoc,thisMartinb,1);
        end
    end
end
    
% Spearman rank correlation coefficient 
matchedSimulatedMartinb = [];
matchedTrueMartinb = [];
corrCoeffLocal = zeros(nLocs,1);
pValueLocal = zeros(nLocs,1);
nMatchupsLocal = zeros(nLocs,1);

for iLoc = 1:nLocs
    martinbTrueLocal = squeeze(martinbTrue(iLoc,:,:));
    martinbSimulatedLocal = squeeze(martinbSimulatedProfile(iLoc,:,:,1));
    validIdxs = ~isnan(martinbTrueLocal) & ~isnan(martinbSimulatedLocal);
    validMartinbTrueLocal = martinbTrueLocal(validIdxs);
    validMartinbSimulatedLocal = martinbSimulatedLocal(validIdxs);

    [corrCoeffLocal(iLoc),pValueLocal(iLoc)] = corr(...
        validMartinbTrueLocal(:),validMartinbSimulatedLocal(:),...
        'Type','Spearman');
    nMatchupsLocal(iLoc) = length(validMartinbTrueLocal);
    matchedSimulatedMartinb = [matchedSimulatedMartinb; validMartinbSimulatedLocal];
    matchedTrueMartinb = [matchedTrueMartinb; validMartinbTrueLocal];

end
[corrCoeffOverall,pValueOverall] = corr(matchedTrueMartinb,matchedSimulatedMartinb,'Type','Spearman');

% R-squared 
rsquared = corrCoeffOverall^2;

% Rest of statistics
[rmse,log_rmse] = calcRootMeanSquaredError(matchedTrueMartinb,matchedSimulatedMartinb); % greek letter: phi
[me,log_me] = calcMeanError(matchedTrueMartinb,matchedSimulatedMartinb); % aka bias, greek letter: delta
[mae,log_mae] = calcMeanAbsoluteError(matchedTrueMartinb,matchedSimulatedMartinb);
mape = calcMeanAbsolutePercentageError(matchedTrueMartinb,matchedSimulatedMartinb);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - PLOT
% -------------------------------------------------------------------------

% Colours
coloursLocs = brewermap(nLocs,'*Set1');

% Axis limits
maxVal = [4, 4, 4, 4, 4, 4];

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.45 0.40],'Color','w') 
haxis = zeros(nLocs,1);

for iSubplot = 1:nLocs

    haxis(iSubplot) = subaxis(2,3,iSubplot,'Spacing',0.028,'Padding',0.028,'Margin',0.05);
    ax(iSubplot).pos = get(haxis(iSubplot),'Position');
    ax(iSubplot).pos(2) = ax(iSubplot).pos(2)+0.030; 
    if (iSubplot == 2 || iSubplot == 5)
        ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.030;
    elseif (iSubplot == 3 || iSubplot == 6)
        ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.060;
    end
    if (iSubplot > 3)
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2)-0.020;
    end
    set(haxis(iSubplot),'Position',ax(iSubplot).pos)

    % Re-order
    switch iSubplot
        case 1
            iLoc = 5; % HOT/ALOHA
        case 2
            iLoc = 4; % BATS/OFP
        case 3
            iLoc = 1; % EqPac
        case 4
            iLoc = 3; % PAP-SO
        case 5
            iLoc = 2; % OSP
        case 6
            iLoc = 6; % HAUSGARTEN
    end

    % Scatter
    x = martinbTrue(iLoc,:,:);
    y = martinbSimulatedProfile(iLoc,:,:,1);
    scatter(haxis(iSubplot),x(:),y(:),30,'o',...
        'MarkerEdgeColor','k','MarkerFaceColor',coloursLocs(iLoc,:),...
        'LineWidth',0.5);
    hold on

    % Set axis limits
    ylim([0 maxVal(iSubplot)])
    xlim([0 maxVal(iSubplot)])
    
    % Plot 1:1 reference line
    hline = refline(haxis(iSubplot),1,0); 
    hline.Color = 'k';
    hline.LineStyle = '--';
    hold off
    
    % Readjust axis limits
    ylim([0 maxVal(iSubplot)])
    xlim([0 maxVal(iSubplot)])

    box on

    % Add statistics
    xt = max(xlim)-0.02*max(xlim); 
    yt = max(ylim)-0.02*max(ylim);
    if (pValueLocal(iLoc) >= 0.05)
        text(xt,yt,... % text position relative to axis
            strcat('{\it r} =',{' '},num2str(corrCoeffLocal(iLoc),'%.2f'),',',...
            {' '},'{\itp} =',{' '},num2str(pValueLocal(iLoc),'%.2f'),',',...
            {' '},'{\itN} =',{' '},num2str(nMatchupsLocal(iLoc),'%.0f')),...
            'FontSize',11,'Horiz','right','Vert','top');
    else
       text(xt,yt,... % text position relative to axis
            strcat('{\it r} =',{' '},num2str(corrCoeffLocal(iLoc),'%.2f'),',',...
            {' '},'{\itp} < 0.05',',',...
            {' '},'{\itN} =',{' '},num2str(nMatchupsLocal(iLoc),'%.0f')),...
            'FontSize',11,'Horiz','right','Vert','top'); 
    end

    if (iSubplot == 5)
        xl = xlabel('True b');
    elseif (iSubplot == 1)
        yl = ylabel('Simulated b');
        yl.Position(2) = yl.Position(2) - 2.7;
        yl.Position(1) = yl.Position(1) - 0.4;
    end

     title(STATION_NAMES(iLoc),'FontSize',14)

end

saveFigure(strcat('samplesize_analysis'))
        
