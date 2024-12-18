
% ======================================================================= %
%                                                                         %
% This script generates Figure 4 in our paper, which compares published   %
% BCP mesopelagic transfer efficiency metrics (Martin's b coefficient, z* % 
% coefficient and Teff) across six ocean sites (HOT/ALOHA, BATS/OFP,      % 
% EqPac, PAP-SO, OSP and HAUSGARTEN). The sites are grouped by ocean      % 
% biome (subtropical, equatorial, subpolar and polar). The script also    %
% produces Figure S6, which shows the monthly evolution of these three    % 
% metrics for each location, comparing data from our sediment             %
% trap/radionuclide compilation with UVP5 data.                           %  
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 26 Nov 2024                                   %
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

% Filename declarations
filenameInputMetricsData           = 'bcpmetrics_all.mat';
filenameInputTimeseriesInformation = 'timeseries_station_information.mat';

% Load the metrics array
load(fullfile('.','data','processed',filenameInputMetricsData),...
    'metricsData','martinbMonthlyTrapAndRad','martinbMonthlyUvp5',...
    'zstarMonthlyTrapAndRad','zstarMonthlyUvp5',...
    'teffMonthlyTrapAndRad','teffMonthlyUvp5')

% Load station information
load(fullfile('.','data','processed',filenameInputTimeseriesInformation),...
    'STATION_NAMES')
nLocs = length(STATION_NAMES);

% Indexes to locations
iE = 1; % EqPac
iO = 2; % OSP
iP = 3; % PAP-SO
iB = 4; % BATS/OFP
iHo = 5; % HOT/ALOHA
iHa = 6; % HAUSGARTEN

% Indexes to metrics
nMetrics = 4;
iMartinb = 1;
iZstar   = 2;
iPeeff   = 3;
iTeff    = 4;

% Indexes to publications from which BCP metrics have been obtained
nPublications = 8 + 2; % +2 to add POC flux compilation and UVP5-derived estimates
iF2002 = 1;  % Francois et al. (2002)
iB2009 = 2;  % Buesseler & Boyd (2009)
iL2011 = 3;  % Lam et al. (2011)
iG2015 = 4;  % Guidi et al. (2015)
iM2016 = 5;  % Mouw et al. (2016)
iH2012 = 6;  % Henson et al. (2012)
iM2015 = 7;  % Marsay et al. (2015)
iW2016 = 8;  % Weber et al. (2016)
iUvp = iW2016 + 1;
iTimeSeries = iUvp + 1;

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - PLOT FIGURE 4
% -------------------------------------------------------------------------

labelMetrics = {'Martin b','z^{*} (m)','T_{eff}'};
labelOceanLocations = {'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'};

nSubplots = nLocs*length(labelMetrics);
nPointsPerLoc = nPublications;
x = 1:nPointsPerLoc;            
myColourPalette = [jet(nPublications-1);[0 0 0]]; % append a row of black at the end

% Array to store group mean
mng = zeros(nPublications,nMetrics,nLocs);

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.95 0.75],'Color','w') 
haxis = zeros(nSubplots,1);

for iSubplot = 1:nSubplots

    % Add an extra row of plots to accommodate the legend
    haxis(iSubplot) = subaxis(4,nLocs,iSubplot,'Spacing',0.020,'Padding',0.024,'Margin',0.02);
    ax(iSubplot).pos = get(haxis(iSubplot),'Position');
    if (iSubplot == 1 || iSubplot == 7 || iSubplot == 13)
        ax(iSubplot).pos(1) = ax(iSubplot).pos(1) + 0.068;
    elseif (iSubplot == 2 || iSubplot == 8 || iSubplot == 14)
        ax(iSubplot).pos(1) = ax(iSubplot).pos(1);
    elseif (iSubplot == 3 || iSubplot == 9 || iSubplot == 15)
        ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.068;
    elseif (iSubplot == 4 || iSubplot == 10 || iSubplot == 16)
        ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.136;
    elseif (iSubplot == 5 || iSubplot == 11 || iSubplot == 17)
        ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.2040;
    elseif (iSubplot == 6 || iSubplot == 12 || iSubplot == 18)
        ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.2720;
    end 
    if (iSubplot >= 1 && iSubplot <= nLocs)
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2); % + 0.02;
    elseif (iSubplot >= nLocs+1 && iSubplot <= 12)
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2); % + 0.01;
    end
    set(haxis(iSubplot),'Position',ax(iSubplot).pos) 

    % Rearrange the order of locations: subtropics - tropics - subpolar
    if (iSubplot == 1 || iSubplot == 7 || iSubplot == 13)
        iLoc = iHo; % HOT/ALOHA
    elseif (iSubplot == 2 || iSubplot == 8 || iSubplot == 14)
        iLoc = iB;  % BATS/OFP
    elseif (iSubplot == 3 || iSubplot == 9 || iSubplot == 15)        
        iLoc = iE;  % EqPac
    elseif (iSubplot == 4 || iSubplot == 10 || iSubplot == 16)         
        iLoc = iP;  % PAP-SO
    elseif (iSubplot == 5 || iSubplot == 11 || iSubplot == 17)
        iLoc = iO;  % OSP
    elseif (iSubplot == 6 || iSubplot == 12 || iSubplot == 18)
        iLoc = iHa; % HAUSGARTEN
    end
    
    % Specify the order of the metrics
    if (iSubplot >= 1 && iSubplot <= 6)
        iMetric = iMartinb; 
    elseif (iSubplot >= 7 && iSubplot <= 12)        
        iMetric = iZstar; 
    elseif (iSubplot >= 13 && iSubplot <= 18)
        iMetric = iTeff; 
    end

    for iRef = 1:nPointsPerLoc 

        % First, plot metrics derived from on-site field observations
        if (iRef < iH2012) 
            
            for iRep = 1:4
                % if the 4th position has 3 values, plot as mean +/- std 
                mny = metricsData(iLoc,iMetric,iRef,1,iRep);
                ypos = metricsData(iLoc,iMetric,iRef,2,iRep) - mny; % length above the data point
                yneg = mny - metricsData(iLoc,iMetric,iRef,3,iRep); % length below the data point
                errorbar(haxis(iSubplot),x(iRef),mny,yneg,ypos,...
                    'o-','Color',myColourPalette(iRef,:),'LineWidth',1.5,...
                    'MarkerSize',4,'CapSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k'); hold on;
            end % iRep
            mng(iRef,iMetric,iLoc) = mean(metricsData(iLoc,iMetric,iRef,1,:),'omitnan');

        % Second, plot metrics derived using some sort of modelling    
        elseif (iRef >= iH2012 && iRef <= iW2016) 

            for iRep = 1:4
                mny = metricsData(iLoc,iMetric,iRef,1,iRep);
                ypos = metricsData(iLoc,iMetric,iRef,2,iRep) - mny; % length above the data point
                yneg = mny - metricsData(iLoc,iMetric,iRef,3,iRep); % length below the data point
                errorbar(haxis(iSubplot),x(iRef),mny,yneg,ypos,...
                    'o-','Color',myColourPalette(iRef,:),'LineWidth',1.5,...
                    'MarkerSize',4,'CapSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k'); hold on;
            end
            mng(iRef,iMetric,iLoc) = mean(metricsData(iLoc,iMetric,iRef,1,:),'omitnan');
        
        % Third, plot UVP5-derived estimates    
        elseif (iRef == nPublications-1) 
        
            mny = metricsData(iLoc,iMetric,iRef,1,1);
            ypos = metricsData(iLoc,iMetric,iRef,2,1) - mny;
            yneg = mny - metricsData(iLoc,iMetric,iRef,3,1);
            errorbar(haxis(iSubplot),x(iRef),mny,yneg,ypos,...
                'o-','Color',myColourPalette(iRef,:),'LineWidth',1.5,...
                'MarkerSize',4,'CapSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k'); hold on;
            mng(iRef,iMetric,iLoc) = mean(metricsData(iLoc,iMetric,iRef,1,:),'omitnan');
        
        % Fourth, plot the metrics calculated from our compilation
        elseif (iRef == nPublications)
                
            mny = metricsData(iLoc,iMetric,iRef,1,1);
            ypos = metricsData(iLoc,iMetric,iRef,2,1) - mny;
            yneg = mny - metricsData(iLoc,iMetric,iRef,3,1);
            errorbar(haxis(iSubplot),x(iRef),mny,yneg,ypos,...
                'o-','Color',myColourPalette(iRef,:),'LineWidth',1.5,...
                'MarkerSize',4,'CapSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k'); hold on;
            mng(iRef,iMetric,iLoc) = mean(metricsData(iLoc,iMetric,iRef,1,:),'omitnan');
            
        end

    end % iRef
    
    % Add group mean and standard deviation
    xvals = cat(2,x(1)-1,x(:)',x(end)+1);
    subplotMean = mean(mng(:,iMetric,iLoc),'omitnan');
    subplotStd = std(mng(:,iMetric,iLoc),'omitnan');
    subplotStdUpp = subplotMean + subplotStd;
    subplotStdLow = subplotMean - subplotStd;
    htmlGray = [128 128 128]/255;
    h19 = plot(haxis(iSubplot),xvals,subplotMean*ones(size(xvals)),...
        'Color',htmlGray,'LineWidth',1.3); hold on;
 
    % Shade the area between the lines
    h20 = fill([xvals, fliplr(xvals)],...
               [subplotMean*ones(size(xvals)),fliplr(subplotStdUpp*ones(size(xvals)))],...
               [0.8, 0.8, 0.8],'FaceAlpha',0.5,'EdgeColor','none');
    h21 = fill([xvals, fliplr(xvals)],...
               [subplotMean*ones(size(xvals)),fliplr(subplotStdLow*ones(size(xvals)))],...
               [0.8, 0.8, 0.8],'FaceAlpha',0.5,'EdgeColor','none');

    uistack(h19,'bottom'); % move to the bottom of the stack of plotted graphical elements
    uistack(h20,'bottom');  
    uistack(h21,'bottom');

    switch iMetric 
        case iMartinb 
            ylim([0 2])
            ytickformat('%.1f')
            yticks([0:0.5:2])
            yticklabels({'0','0.5','1.0','1.5','2.0'});
        case iZstar 
            ylim([10 1250]) 
            yticks([0:250:1250])
            yticklabels({'0','250','500','750','1000','1250'});
        case iTeff 
            ylim([0 0.80]) 
            ytickformat('%.1f')
            yticks([0:0.20:0.80])
            yticklabels({'0','0.20','0.40','0.60','0.80'});
    end
     
    if (iSubplot == 1)
        ylabel(labelMetrics(1),'FontSize',12,'Units','normalized','Position',[-0.3, 0.5, 0])
    elseif (iSubplot == 7)        
        ylabel(labelMetrics(2),'FontSize',12,'Units','normalized','Position',[-0.3, 0.5, 0])
    elseif (iSubplot == 13)
        ylabel(labelMetrics(3),'FontSize',12,'Units','normalized','Position',[-0.3, 0.5, 0])
    else
        set(gca,'yticklabel',[])
    end
    
    % Add grid
    ah = gca;
    ah.XAxis.FontSize = 7;
    ah.XGrid = 'off';

    % xlabels
    xlim([-0.4 nPublications+1.4])
    xticks(1:nPublications)
    xticklabels({'Fr2002',...
                'Bu2009',...
                'La2011',...
                'Gu2015',...
                'Mo2016',...
                'He2012',...                  
                'Ma2015',...
                'We2016',...
                'UVP5',...
                'T&R'});
    xtickangle(90)

    % Tune box and grid
    ah.Box = 'off';
    if (iSubplot == 1 || iSubplot == 7 || iSubplot == 13)
        ah.YColor = 'k';
    else
        ah.YColor = [0.7, 0.7, 0.7, 0.3]; % RGBA: [1, 1, 1] for white, 0 for 100% transparency
%         ah.YColor = 'w'; 
    end
    ah.YGrid = 'on';
    ah.GridColor = 'k';  

    if (iSubplot >= 1 && iSubplot <= nLocs)
        title(labelOceanLocations(iSubplot),'FontSize',12,'Units','normalized','Position',[0.5, 1.05, 0])
    end
    
    % Add legend
    if (iSubplot == nSubplots)
 
        for iRef = 1:nPublications
            qw{iRef} = scatter(nan,50,'o','MarkerEdgeColor',myColourPalette(iRef,:),...
                'MarkerFaceColor',myColourPalette(iRef,:),'LineWidth',1.2);
        end
        lg = legend([qw{:}], {'Francois et al. (2002)',...
                              'Buesseler & Boyd (2009)',...
                              'Lam et al. (2011)',...
                              'Guidi et al. (2015)',...
                              'Mouw et al. (2016b)',...
                              'Henson et al. (2012)',...                  
                              'Marsay et al. (2015)',...
                              'Weber et al. (2016)',...
                              'UVP5 compilation (Kiko et al., 2022)',...
                              'Trap & radionuclide compilation (this study)'},...
                              'NumColumns',2);
        lg.Position(1) = 0.26; lg.Position(2) = 0.12;
        lg.Orientation = 'vertical';
        lg.FontSize = 12; 
        lg.ItemTokenSize = [24,1];
        set(lg,'Box','off')   
        
    end     
       
end % iSubplot

saveFigure('uncertainty_bcpmetrics')

clear ax

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - PLOT FIGURE S6
% -------------------------------------------------------------------------

coloursPlots = [1, 0, 0; 0, 0, 0]; % red for UVP, black for sediment trap/radionuclide

maxValMetric = [3, 3, 4, 3, 3, 3;...
                2100, 2100, 2100, 2100, 2100, 2100;...
                1.2, 1.2, 1.2, 1.2, 1.2, 3.5];

% Axis limits
monthLabels = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
metricTags = {'martinb','zstar','teff'};

for iFigure = 1:3
    
    switch iFigure
        case 1
            iMetric = iMartinb;
            titleStr = {'Martin b'};
            metricTrapAndRad_perm = permute(martinbMonthlyTrapAndRad, [2 1 3]); % 12 months x nLocs x 3 (mean, upp, low)
            metricUvp5_perm       = permute(martinbMonthlyUvp5, [2 1 3]);       % 12 months x nLocs x 3 (mean, upp, low)
        case 2
            iMetric = iZstar;
            titleStr = {'z*'};
            metricTrapAndRad_perm = permute(zstarMonthlyTrapAndRad, [2 1 3]); % 12 months x nLocs x 3 (mean, upp, low)
            metricUvp5_perm       = permute(zstarMonthlyUvp5, [2 1 3]);       % 12 months x nLocs x 3 (mean, upp, low)
        case 3
            iMetric = iTeff;
            titleStr = {'Teff'};
            metricTrapAndRad_perm = permute(teffMonthlyTrapAndRad, [2 1 3]); % 12 months x nLocs x 3 (mean, upp, low)
            metricUvp5_perm       = permute(teffMonthlyUvp5, [2 1 3]);       % 12 months x nLocs x 3 (mean, upp, low)
    end
    
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.55 0.45],'Color','w') 
    haxis = zeros(nLocs,1);

    for iSubplot = 1:nLocs

        haxis(iSubplot) = subaxis(2,3,iSubplot,'Spacing',0.028,'Padding',0.028,'Margin',0.10);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');
        ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.035;
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2)+0.030; 
        if (iSubplot == 2 || iSubplot == 5)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.035;
        elseif (iSubplot == 3 || iSubplot == 6)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.070;
        end
        if (iSubplot > 3)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2)-0.040;
        end
        set(haxis(iSubplot),'Position',ax(iSubplot).pos)

        % Re-order
        switch iSubplot
            case 1
                iLoc = iHo; % HOT/ALOHA
            case 2
                iLoc = iB; % BATS/OFP
            case 3
                iLoc = iE; % EqPac
            case 4
                iLoc = iP; % PAP-SO
            case 5
                iLoc = iO; % OSP
            case 6
                iLoc = iHa; % HAUSGARTEN
        end

        vals1 = metricUvp5_perm(:,iLoc,1);
        vals2 = metricTrapAndRad_perm(:,iLoc,1);

        pos1 = metricUvp5_perm(:,iLoc,2) - vals1;
        pos2 = metricTrapAndRad_perm(:,iLoc,2) - vals2;

        neg1 = vals1 - metricUvp5_perm(:,iLoc,3);
        neg2 = vals2 - metricTrapAndRad_perm(:,iLoc,3);

        % Grouped bar plot
        hbar = bar(haxis(iSubplot),[vals1, vals2],'grouped','BarWidth',1,'FaceColor','flat'); 
        hold on

        % Colour bars
        for k = 1:numel(hbar)
            hbar(k).CData = coloursPlots(k,:);
        end
        hold on

        % Add error bars
        for k = 1:numel(hbar)                                                      
            xtips = hbar(k).XEndPoints;
            ytips = hbar(k).YEndPoints;
            if k == 1
                errorbar(haxis(iSubplot),xtips,ytips,neg1,pos1,'.k','MarkerSize',0.2,'CapSize', 2,'HandleVisibility','off')
            elseif k == 2
                errorbar(haxis(iSubplot),xtips,ytips,neg2,pos2,'.k','MarkerSize',0.2,'CapSize', 2,'HandleVisibility','off')
            end
            hold on
        end
        hold on

        % Add annual mean from UVP5 compilation
        yline(metricsData(iLoc,iMetric,9,1,1),'-','Color',coloursPlots(1,:),'LineWidth',1);
        hold off

        % Add annual mean from sediment trap and radionuclide compilation
        yline(metricsData(iLoc,iMetric,10,1,1),'-','Color',coloursPlots(2,:),'LineWidth',1);
        hold on

        % Set yaxis limits and label
        ylim([0 maxValMetric(iFigure,iSubplot)]);
        if (iFigure == 1 || iFigure ==3)
            ytickformat('%.1f')
        else
            ytickformat('%.0f')
        end

        % Set xaxis labels
        set(gca,'xticklabel',monthLabels);
        xtickangle(90);

        % Grid
        axgrid = gca;
        axgrid.YGrid = 'on';  % Enable only horizontal grid lines
        axgrid.XGrid = 'off'; % Disable vertical grid lines

        % Title
        title(STATION_NAMES(iLoc),'FontSize',12)
        
        % Legend
        if (iSubplot == nLocs)
            lg = legend('UVP5','T&R','Annual mean UVP5','Annual mean T&R','Location','eastoutside');
            lg.Position(1) = 0.80; lg.Position(2) = 0.78;
            lg.Orientation = 'vertical';
            lg.FontSize = 12; 
            lg.ItemTokenSize = [20,5];
            set(lg,'Box','off')   
        end

    end % iSubplot

    % Give common title to the figure
    a = axes;
    t = title(titleStr,'FontSize',16);
    % Specify visibility of the current axis as 'off'
    a.Visible = 'off';
    % Specify visibility of Title, XLabel, and YLabel as 'on'
    t.Visible = 'on';
    t.Position(1) = t.Position(1);
    t.Position(2) = t.Position(2) + 0.040;

    saveFigure(strcat('barplot_',metricTags{iFigure}))

end % iFigure
