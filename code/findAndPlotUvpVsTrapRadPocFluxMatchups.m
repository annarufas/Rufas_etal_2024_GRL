
% ======================================================================= %
%                                                                         %
% This script finds matchups between the POC flux data from the trap and  %
% radionuclide compilation and the POC flux estimated from the UVP5       % 
% dataset. The script is structured into 4 sections:                      %
%   Section 1 - Presets.                                                  %           
%   Section 2 - Find matchups between the UVP5-derived estimates and the  %
%               trap and radionuclide measurements.                       %
%   Section 3 - Compute correlation coefficients.                         %
%   Section 4 – Plot matchups figure (Figure 3).                          %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 14 Oct 2024                                   %
%                                                                         %
% ======================================================================= %

close all; clear all; clc
addpath(genpath('./data/raw/'));
addpath(genpath('./data/processed/'));
addpath(genpath('./code/'));
addpath(genpath('./resources/external/'));
addpath(genpath('./resources/internal/'));

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

filenameInputUvpProcessedDataset45sc = 'pocflux_bisson_45sc.mat';
filenameInputPocFluxCompilation      = 'pocflux_compilation.mat';
filenameInputTimeseriesInformation   = 'timeseries_station_information.mat';

% Load the UVP5 data
load(fullfile('.','data','processed','UVP5',filenameInputUvpProcessedDataset45sc),'UVP_TABLE')

% Load the trap and radionuclide compilation
load(fullfile('.','data','processed',filenameInputPocFluxCompilation),'TRAPRAD_TABLE')

% Load station information
load(fullfile('.','data','processed',filenameInputTimeseriesInformation),'STATION_NAMES')
nLocs = length(STATION_NAMES);
SUFFIX_ECOTAXA_FOLDER_NAME = {'EqPac','OSP','PAPSO','BATSOFP','HOTALOHA','HAUSGARTEN'};

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - FIND MATCHUPS BETWEEN THE UVP5-DERIVED VALUES AND THE TRAP
% AND RADIONUCLIDE MEASUREMENTS
% -------------------------------------------------------------------------

% Define the tolerance values
depthTolerance = 5; % 5 metres
timeTolerance = 2;  % 2 days

% Initialise structure to hold data for plotting
MATCHUP_TABLE = struct();

% Initialise counters
nMatchups = 0;
nAvailableUvpSamplesForMatchup = 0;
nMatchedLocs = 0;

% Initialise cell array that contains names of locations with matchups
tableNamesMatchedLocs = {};

for iLoc = 1:nLocs
    
    thisLocationName = SUFFIX_ECOTAXA_FOLDER_NAME{iLoc};

    % Extract local data
    currLocUvp = UVP_TABLE(strcmp(UVP_TABLE.Location,SUFFIX_ECOTAXA_FOLDER_NAME{iLoc}),:);
    currLocTraprad = TRAPRAD_TABLE(TRAPRAD_TABLE.tag == STATION_NAMES{iLoc},:);
    
    % Calculate net error for traprad data
    currLocTraprad.netError = sqrt(currLocTraprad.randerr_POC_mg_m2_d.^2 +... 
        (currLocTraprad.syserr_POC_mmol_m2_d.*12.011).^2);

    % The UVP data has several casts in a day. Let's group data taken in 
    % same depths and day and calculate the average and propagate error 
    % before attempting matchups
    currLocUvp.DayOnly = dateshift(currLocUvp.yyyy_mm_ddHh_mm, 'start', 'day'); % convert datetime to just date (ignores time)
    [G, zz, tt] = findgroups(currLocUvp.Depth_m_, currLocUvp.DayOnly); % find groups based on depth and time
    avgPocFluxUvp = splitapply(@mean, currLocUvp.modelledPocFlux, G); % calculate the mean POC flux for each group
    
    % Loop through each unique group to handle error
    nGroups = max(G); % no. unique groups
    errPocFluxUvp = zeros(nGroups, 1); 
    for i = 1:nGroups
        vals = currLocUvp.modelledPocFlux(G == i);
        errs = currLocUvp.fluxErr(G == i);
        % Apply the worstcase function
        [x_LB,x_UB,f_LB,f_MID,f_UB,minus_percent,plus_percent] = ... 
            worstcase(@(x) mean(x),vals,errs);
        errPocFluxUvp(i) = f_UB - f_MID; 
    end

    % Create a new table with the grouped results
    currLocUvpAvgPocFlux = table(zz, tt, avgPocFluxUvp, errPocFluxUvp); 
    currLocUvpAvgPocFlux.tt = datestr(currLocUvpAvgPocFlux.tt, 'yyyy-mm-dd'); 
    nAvailableUvpSamplesForMatchup = nAvailableUvpSamplesForMatchup + height(currLocUvpAvgPocFlux);
    
    % Convert date in currLocTraprad to match date format in currLocUvpAvgPocFlux
    dateObj = datetime(currLocTraprad.midDate, 'InputFormat', 'dd-MMM-yyyy');
    currLocTraprad.midDate = datestr(dateObj, 'yyyy-mm-dd');

    % Convert dates to numeric formats
    uvp_dates_numeric = datenum(currLocUvpAvgPocFlux.tt);
    traprad_dates_numeric = datenum(currLocTraprad.midDate);

    % Loop through each row in currLocTraprad and find matching rows in currLocUvpAvgPocFlux
    matchedIndices = [];
    
    for i = 1:height(currLocTraprad)

        % Calculate the absolute difference in depth and date
        depthDiff = abs(currLocUvpAvgPocFlux.zz - currLocTraprad.depth(i));
        dateDiff = abs(uvp_dates_numeric - traprad_dates_numeric(i));

        % Find matches where both depth and date differences are within tolerance
        isMatch = (depthDiff <= depthTolerance) & (dateDiff <= timeTolerance);

        % Get indices of matches
        matchingIndices = find(isMatch);
        
        % To handle cases where more than one match is found, you can select 
        % the match that is closest to both the depth and time by calculating 
        % a "distance" measure (a combined score based on depth and time 
        % differences). This way, you can choose the match with the smallest 
        % combined difference.
        if ~isempty(matchingIndices)
            if length(matchingIndices) > 1
                combinedDiff = depthDiff(matchingIndices) + dateDiff(matchingIndices);
                minVal = min(combinedDiff);
                bestMatchIdx = matchingIndices(combinedDiff == minVal);
                for j = 1:length(bestMatchIdx)
                    matchedIndices = [matchedIndices; i, bestMatchIdx(j)]; 
                end
            else
                matchedIndices = [matchedIndices; i, matchingIndices];
            end
        end

    end
    
    if ~isempty(matchedIndices)
        
        % Extract the matching rows from both tables based on the matched indices
        matchedTrapRadPocFlux    = currLocTraprad.POC_mg_m2_d(matchedIndices(:,1));
        matchedTrapRadPocFluxErr = currLocTraprad.netError(matchedIndices(:,1));
        matchedUvpPocFlux        = currLocUvpAvgPocFlux.avgPocFluxUvp(matchedIndices(:,2));
        matchedUvpPocFluxErr     = currLocUvpAvgPocFlux.errPocFluxUvp(matchedIndices(:,2));

        % Add to counter
        nMatchups = nMatchups + height(matchedTrapRadPocFlux);

        % Store
        MATCHUP_TABLE.(thisLocationName).traprad.avg = matchedTrapRadPocFlux;
        MATCHUP_TABLE.(thisLocationName).traprad.err = matchedTrapRadPocFluxErr;
        MATCHUP_TABLE.(thisLocationName).uvp.avg = matchedUvpPocFlux;
        MATCHUP_TABLE.(thisLocationName).uvp.err = matchedUvpPocFluxErr;

        nMatchedLocs = nMatchedLocs + 1;
        tableNamesMatchedLocs = [tableNamesMatchedLocs; thisLocationName]; 
    
    end
    
end

fprintf('%d matchups were found between measured and estimated values of POC flux.\n', nMatchups)
fprintf('%d UVP5 samples were available.\n', nAvailableUvpSamplesForMatchup)

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - COMPUTE CORRELATION COEFFICIENTS
% -------------------------------------------------------------------------

% Spearman rank correlation coefficient 
matchedUvpPocFlux = [];
matchedTrapRadPocFlux = [];
corrCoeffLocal = zeros(nMatchedLocs,1);
pValueLocal = zeros(nMatchedLocs,1);
nMatchupsLocal = zeros(nMatchedLocs,1);

iMatchedLoc = 0;
for iLoc = 1:nLocs
    thisLocationName = SUFFIX_ECOTAXA_FOLDER_NAME{iLoc};
    if isfield(MATCHUP_TABLE, thisLocationName)
        iMatchedLoc = iMatchedLoc + 1;
        [corrCoeffLocal(iMatchedLoc),pValueLocal(iMatchedLoc)] = corr(...
            MATCHUP_TABLE.(thisLocationName).traprad.avg,...
            MATCHUP_TABLE.(thisLocationName).uvp.avg,...
            'Type','Spearman');
        nMatchupsLocal(iMatchedLoc) = length(MATCHUP_TABLE.(thisLocationName).uvp.avg);
        matchedUvpPocFlux = [matchedUvpPocFlux; MATCHUP_TABLE.(thisLocationName).uvp.avg];
        matchedTrapRadPocFlux = [matchedTrapRadPocFlux; MATCHUP_TABLE.(thisLocationName).traprad.avg];
    end 
end
[corrCoeffOverall,pValueOverall] = corr(matchedTrapRadPocFlux,matchedUvpPocFlux,'Type','Spearman');

% R-squared 
rsquared = corrCoeffOverall^2;

% Rest of statistics
[rmse,log_rmse] = calcRootMeanSquaredError(matchedTrapRadPocFlux,matchedUvpPocFlux); % greek letter: phi
[me,log_me] = calcMeanError(matchedTrapRadPocFlux,matchedUvpPocFlux); % aka bias, greek letter: delta
[mae,log_mae] = calcMeanAbsoluteError(matchedTrapRadPocFlux,matchedUvpPocFlux);
mape = calcMeanAbsolutePercentageError(matchedTrapRadPocFlux,matchedUvpPocFlux);
       
% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - PLOT MATCHUPS FIGURE (FIGURE 3)
% -------------------------------------------------------------------------

% .........................................................................

% Plot all

coloursLocs = brewermap(nMatchedLocs,'*Set1');

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.28 0.36],'Color','w')
axh = axes('Position', [0.11 0.11 0.62 0.72]); % axes to make space for the legend on the right-hand-side

% Plot error bars first
iMatchedLoc = 0;
for iLoc = 1:nLocs
    thisLocationName = SUFFIX_ECOTAXA_FOLDER_NAME{iLoc};
    
    if isfield(MATCHUP_TABLE, thisLocationName)
        
        x = MATCHUP_TABLE.(thisLocationName).traprad.avg;
        y = MATCHUP_TABLE.(thisLocationName).uvp.avg;
        xerr = MATCHUP_TABLE.(thisLocationName).traprad.err;
        yerr = MATCHUP_TABLE.(thisLocationName).uvp.err;
        iMatchedLoc = iMatchedLoc + 1;

        scatter(axh,NaN,NaN,60,'o','MarkerEdgeColor','k','MarkerFaceColor',coloursLocs(iMatchedLoc,:),...
            'LineWidth',0.5,'HandleVisibility','off');
        hold on
        
        eb(1) = errorbar(x,y,xerr,'horizontal',...
            'LineStyle','none','HandleVisibility','off');
        eb(2) = errorbar(x,y,yerr,'vertical',...
            'LineStyle','none','HandleVisibility','off');
        eb(1).CapSize = 0;
        eb(2).CapSize = 0;
        set(eb,'Color',coloursLocs(iMatchedLoc,:),'LineWidth',1)

        % Set transparency level (0:1)
        alpha = 0.3;  
        set([eb(1).Bar,eb(1).Line],'ColorType','truecoloralpha',...
            'ColorData',[eb(1).Line.ColorData(1:3); 255*alpha])
        set([eb(2).Bar,eb(2).Line],'ColorType','truecoloralpha',...
            'ColorData',[eb(2).Line.ColorData(1:3); 255*alpha])
        hold on

    end
end
hold on

% Plot scatters on top
iMatchedLoc = 0;
for iLoc = 1:nLocs  
    thisLocationName = SUFFIX_ECOTAXA_FOLDER_NAME{iLoc};
    
    if isfield(MATCHUP_TABLE, thisLocationName)  
        x = MATCHUP_TABLE.(thisLocationName).traprad.avg;
        y = MATCHUP_TABLE.(thisLocationName).uvp.avg;
        iMatchedLoc = iMatchedLoc + 1;

        scatter(axh,x,y,60,'o','MarkerEdgeColor','k',...
            'MarkerFaceColor',coloursLocs(iMatchedLoc,:),'LineWidth',0.5);
        hold on
    end
    
    % Plot 1:1 reference line
    if (iLoc == nLocs) 
        hline = refline(axh,1,0); 
        hline.Color = 'k';
        hline.LineStyle = '--';
    end
    hold on

end
hold on
box on

ylim([0 900])
xlim([0 220]) 

% Add statistics
xt = max(xlim)-0.02*max(xlim); 
yt = max(ylim)-0.02*max(ylim);
if (pValueOverall >= 0.01)
    text(xt,yt,... % text position relative to axis
        strcat('{\it r} =',{' '},num2str(corrCoeffOverall,'%.2f'),',',...
        {' '},'{\itp} =',{' '},num2str(pValueOverall,'%.2f'),',',...
        {' '},'{\itN} =',{' '},num2str(nMatchups,'%.0f')),...
        'FontSize',12,'Horiz','right','Vert','top','FontWeight','bold');
else
   text(xt,yt,... % text position relative to axis
        strcat('{\it r} =',{' '},num2str(corrCoeffOverall,'%.2f'),',',...
        {' '},'{\itp} < 0.01',...
        {' '},'{\itN} =',{' '},num2str(nMatchups,'%.0f')),...
        'FontSize',12,'Horiz','right','Vert','top','FontWeight','bold'); 
end

% Add legend
lg = legend(axh,{'OSP','PAP-SO','BATS/OFP','HOT/ALOHA','1:1 line'});  
lg.Position(1) = 0.74; lg.Position(2) = 0.63;
lg.ItemTokenSize = [15,1];
lg.FontSize = 11;
set(lg,'Box','off')

yl = ylabel('Estimated (UVP5)');
xl = xlabel('Measured (sediment traps & radionuclides)');
yl.Position(1) = yl.Position(1) - 3;
xl.Position(2) = xl.Position(2) - 10;

set(axh,'FontSize',12)

saveFigure('matchups_estimated_vs_measured_45sizeclasses_all')

% .........................................................................

% Plot by station

% Adjust these accordingly (matches the re-order criterion)
xMax = [40, 130, 142, 93];
yMax = [78, 350, 939, 45];

coloursLocs = brewermap(nMatchedLocs,'*Set1');

% Adjust accordingly (matches the re-order criterion)
titlesMatchedLocs = {'HOT/ALOHA','BATS/OFP','PAP-SO','OSP'};

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.31 0.46],'Color','w')
haxis = zeros(nMatchedLocs,1);

for iLoc = 1:nMatchedLocs

    % Re-order
    switch iLoc
        case 1
            iMatchedLoc = 4; % HOT/ALOHA
        case 2
            iMatchedLoc = 3; % BATS/OFP
        case 3
            iMatchedLoc = 2; % PAP-SO
        case 4
            iMatchedLoc = 1; % OSP
    end 
    
    thisLocationName = tableNamesMatchedLocs{iMatchedLoc};

    haxis(iLoc) = subaxis(2,2,iLoc,'Spacing',0.02,'Padding',0.04,'Margin',0.05);
    ax(iLoc).pos = get(haxis(iLoc),'Position');
    if (iLoc == 3 || iLoc == 4)
        ax(iLoc).pos(2) = ax(iLoc).pos(2)+0.010; 
    end
    set(haxis(iLoc),'Position',ax(iLoc).pos)

    x = MATCHUP_TABLE.(thisLocationName).traprad.avg;
    y = MATCHUP_TABLE.(thisLocationName).uvp.avg;
    xerr = MATCHUP_TABLE.(thisLocationName).traprad.err;
    yerr = MATCHUP_TABLE.(thisLocationName).uvp.err;

    % Plot error bars first
    scatter(haxis(iLoc),NaN,NaN,60,'o','MarkerEdgeColor','k','MarkerFaceColor',coloursLocs(iLoc,:),...
        'LineWidth',0.5,'HandleVisibility','off');
    hold on
    eb(1) = errorbar(x,y,xerr,'horizontal',...
        'LineStyle','none','HandleVisibility','off');
    eb(2) = errorbar(x,y,yerr,'vertical',...
        'LineStyle','none','HandleVisibility','off');
    eb(1).CapSize = 0;
    eb(2).CapSize = 0;
    set(eb,'Color',coloursLocs(iLoc,:),'LineWidth',1)

    % Set transparency level for error bars (0:1)
    alpha = 0.3;  
    set([eb(1).Bar,eb(1).Line],'ColorType','truecoloralpha',...
        'ColorData',[eb(1).Line.ColorData(1:3); 255*alpha])
    set([eb(2).Bar,eb(2).Line],'ColorType','truecoloralpha',...
        'ColorData',[eb(2).Line.ColorData(1:3); 255*alpha])
    hold on

    % Plot scatters on top
    scatter(haxis(iLoc),x,y,60,'o','MarkerEdgeColor','k',...
        'MarkerFaceColor',coloursLocs(iLoc,:),'LineWidth',0.5);
    hold on

    % Define axis limits
    ylim([0 yMax(iLoc)])
    xlim([0 xMax(iLoc)])

    % Plot 1:1 reference line
    hline = refline(haxis(iLoc),1,0); 
    hline.Color = 'k';
    hline.LineStyle = '--';
    hold off
    
    % Readjust axis limits
    ylim([0 yMax(iLoc)])
    xlim([0 xMax(iLoc)])

    box on

    % Add statistics
    xt = max(xlim)-0.02*max(xlim); 
    yt = max(ylim)-0.02*max(ylim);
    if (pValueLocal(iMatchedLoc) >= 0.01)
        text(xt,yt,... % text position relative to axis
            strcat('{\it r} =',{' '},num2str(corrCoeffLocal(iMatchedLoc),'%.2f'),',',...
            {' '},'{\itp} =',{' '},num2str(pValueLocal(iMatchedLoc),'%.2f'),',',...
            {' '},'{\itN} =',{' '},num2str(nMatchupsLocal(iMatchedLoc),'%.0f')),...
            'FontSize',11,'Horiz','right','Vert','top');
    else
       text(xt,yt,... % text position relative to axis
            strcat('{\it r} =',{' '},num2str(corrCoeffLocal(iMatchedLoc),'%.2f'),',',...
            {' '},'{\itp} < 0.01',',',...
            {' '},'{\itN} =',{' '},num2str(nMatchupsLocal(iMatchedLoc),'%.0f')),...
            'FontSize',11,'Horiz','right','Vert','top'); 
    end

    if (iLoc == nMatchedLocs)
        yl = ylabel('Estimated (UVP5)');
        xl = xlabel('Measured (sediment traps & radionuclides)');
        yl.Position(1) = yl.Position(1) - 123;
        yl.Position(2) = yl.Position(2) + 24;
        xl.Position(2) = xl.Position(2) - 3;
        xl.Position(1) = xl.Position(1) - 53;
    end

    title(titlesMatchedLocs{iLoc})
   
end % iLoc

saveFigure('matchups_estimated_vs_measured_45sizeclasses_bystat')
