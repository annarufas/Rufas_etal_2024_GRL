
% ======================================================================= %
%                                                                         %
% This script produces all the figures in our paper related to our POC    %
% flux compilation of sediment trap and radionuclide data: Figure 2 in    % 
% the main text, and Figures S1 and S2 in the Supplementary. The figures  %
% are finished on PPT. The script has 5 sections:                         %
%   Section 1 - Presets.                                                  %
%   Section 2 - Plot Figure 2 (POC flux data by station and month).       %
%   Section 3 - Plot Figure S1 (no. entries by month, location and depth  %
%               horizon).                                                 %
%   Section 4 - Plot percentage uncertainty by month, location and depth  %        
%               horizon.                                                  %
%   Section 5 - Plot Figure S2 (monthly fluxes and their error by         %
%               location and depth horizon).                              %
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

filenameInputMonthlyPocFlux = 'pocflux_compilation_monthly.mat';
filenameInputTimeseriesInformation = 'timeseries_station_information.mat';

load(fullfile('.','data','processed',filenameInputMonthlyPocFlux),...
    'classicRawProfileValues','classicRawProfileDepths','classicRawProfileDataType',...
    'classicRawDhValues_cell','classicRawDhDepths_cell','classicRawDhDataType_cell',...
    'classicMonthlyDhAvg','classicMonthlyDhN','classicMonthlyDhErrTot')

load(fullfile('.','data','processed',filenameInputTimeseriesInformation),...
    'LOC_DEPTH_HORIZONS','STATION_NAMES','STATION_TAGS','NUM_LOCS','NUM_TARGET_DEPTHS')

% Parameters
MOLAR_MASS_CARBON = 12.011; % g mol-1
MAX_NUM_VALUES_PER_MONTH = 1000;

monthLabel = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - PLOT FIGURE 2 (POC FLUX DATA BY STATION AND MONTH)
% -------------------------------------------------------------------------

% localMaxPocFlux = zeros(NUM_LOCS,12);
% for iLoc = 1:NUM_LOCS
%     for iMonth = 1:12
%         allmyvals = squeeze(classicRawProfileValues(:,iMonth,iLoc));
%         if (~isempty(allmyvals))
%             localMaxPocFlux(iLoc,iMonth) = MOLAR_MASS_CARBON.*max(allmyvals); % mmol C m-2 d-1 --> mg C m-2 d-1
%         end
%     end
%     disp(max(localMaxPocFlux(iLoc,:)))
% end

mycolours = parula(3);

for iLoc = 1:NUM_LOCS

    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.25 0.75],'Color','w')
    haxis = zeros(12,1);
    
    for iMonth = 1:12

        haxis(iMonth) = subaxis(4,3,iMonth,'Spacing',0.055,'Padding',0,'Margin',0.13);
        ax(iMonth).pos = get(haxis(iMonth),'Position');

        %%%%%%%%%%%%%% ALL DATA %%%%%%%%%%%%%%
        
        yall = squeeze(classicRawProfileValues(:,iMonth,iLoc)).*MOLAR_MASS_CARBON;
        xall = squeeze(classicRawProfileDepths(:,iMonth,iLoc));
        gall = squeeze(classicRawProfileDataType(:,iMonth,iLoc));

        %%%%%%%%%%%%%% DH DATA %%%%%%%%%%%%%%

        currMonthDhData     = zeros(NUM_TARGET_DEPTHS,MAX_NUM_VALUES_PER_MONTH);
        currMonthDhDataType = cell(NUM_TARGET_DEPTHS,MAX_NUM_VALUES_PER_MONTH);
        currMonthDhDepths   = zeros(NUM_TARGET_DEPTHS,MAX_NUM_VALUES_PER_MONTH);
        currMonthDhN        = zeros(NUM_TARGET_DEPTHS,1);

        for iDh = 1:NUM_TARGET_DEPTHS
            allmyvals      = classicRawDhValues_cell{iDh,iMonth,iLoc};
            allmydepths    = classicRawDhDepths_cell{iDh,iMonth,iLoc};
            allmydatatypes = classicRawDhDataType_cell{iDh,iMonth,iLoc};
            if (~isempty(allmyvals))
                currMonthDhN(iDh) = length(str2num(allmyvals));
                currMonthDhData(iDh,1:currMonthDhN(iDh)) = MOLAR_MASS_CARBON.*str2num(allmyvals); % mmol C m-2 d-1 --> mg C m-2 d-1
                thetypes = textscan(allmydatatypes,'%s');
                for iDataPoint = 1:currMonthDhN(iDh)
                    currMonthDhDataType{iDh,iDataPoint} = thetypes{1}{iDataPoint};
                end
                currMonthDhDepths(iDh,1:currMonthDhN(iDh)) = str2num(allmydepths);
            end
        end

        [rowIdxs,colIdxs,y] = find(currMonthDhData);
        nDataPointsInDh = length(y);
        x = zeros(nDataPointsInDh,1);
        g = cell(nDataPointsInDh,1);
        for iDataPoint = 1:nDataPointsInDh
            x(iDataPoint) = currMonthDhDepths(rowIdxs(iDataPoint),colIdxs(iDataPoint));
            g(iDataPoint) = currMonthDhDataType(rowIdxs(iDataPoint),colIdxs(iDataPoint));
            if (strcmp(g{iDataPoint},'trap'))
                if (x(iDataPoint) > LOC_DEPTH_HORIZONS(iLoc,1,1)... 
                        && x(iDataPoint) < LOC_DEPTH_HORIZONS(iLoc,2,1))
                    g(iDataPoint) = {'zeu'};
                elseif (x(iDataPoint) > LOC_DEPTH_HORIZONS(iLoc,1,2)... 
                        && x(iDataPoint) < LOC_DEPTH_HORIZONS(iLoc,2,2))
                    g(iDataPoint) = {'zmeso'};
                elseif (x(iDataPoint) > LOC_DEPTH_HORIZONS(iLoc,1,3)... 
                        && x(iDataPoint) < LOC_DEPTH_HORIZONS(iLoc,2,3))
                    g(iDataPoint) = {'zbathy'};
                end
            end
        end
        
        % .................................................................

        % We have 4 possible combinations. Mask the entries of interest.
        % If there are no entries for the specific combination 
        % (e.g., zeu-radioisotope), plot NaN, so that the plot generates 
        % an entry for that combination
            
        % All observations

        % Sediment trap
        mask = strcmp(gall,'trap');
        if (sum(mask) == 0)
            d01 = plot(NaN, NaN, 'o', 'MarkerEdgeColor', [0.85 0.85 0.85],...
                'MarkerFaceColor', [0.8 0.8 0.8], 'LineWidth', 1.5,... 
                'HandleVisibility', 'off');
        else
            d01 = plot(yall(mask), xall(mask), 'o', 'MarkerEdgeColor', [0.85 0.85 0.85],...
                'MarkerFaceColor', [0.8 0.8 0.8], 'LineWidth', 1.5,... 
                'HandleVisibility', 'off');
        end
        max01 = max(yall(mask));
        hold on
        % Radionuclide
        mask = strcmp(gall,'radionuclide');
        if (sum(mask) == 0)
            d02 = plot(NaN, NaN, '+', 'MarkerEdgeColor', [0.85 0.85 0.85],...
                'MarkerFaceColor', [0.8 0.8 0.8], 'LineWidth', 1.5,...
                'HandleVisibility', 'off');
        else
            d02 = plot(yall(mask), xall(mask), '+', 'MarkerEdgeColor', [0.85 0.85 0.85],...
                'MarkerFaceColor', [0.8 0.8 0.8], 'LineWidth', 1.5,... 
                'HandleVisibility', 'off');
        end
        max02 = max(yall(mask));
        hold on

        % Observations by depth horizon (zeu, zmeso and zbathy)

        % zeu - sediment traps
        mask = strcmp(g,'zeu');
        if (sum(mask) == 0)
            d1 = plot(NaN, NaN, 'o','MarkerEdgeColor', 'k',...
                'MarkerFaceColor', mycolours(3,:), 'Linewidth', 0.5,...
                'DisplayName', 'trap, z_{eu}');
        else
            d1 = plot(y(mask), x(mask), 'o', 'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', mycolours(3,:), 'LineWidth', 0.5,...
                'DisplayName', 'trap, z_{eu}');
        end
        max1 = max(y(mask));
        hold on
        % zeu - radionuclides
        mask = strcmp(g,'radionuclide');
        if (sum(mask) == 0)
            d2 = plot(NaN, NaN, '+k', 'LineWidth', 1.5,... 
                'DisplayName', 'radionuclide, z_{eu}');
        else
            d2 = plot(y(mask), x(mask), '+k', 'LineWidth', 1.5,...
                'DisplayName', 'radionuclide, z_{eu}');
        end
        max2 = max(y(mask));
        hold on
        % zmeso - sediment traps
        mask = strcmp(g,'zmeso');
        if (sum(mask) == 0)
            d3 = plot(NaN, NaN, 'o', 'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', mycolours(2,:), 'LineWidth', 0.5,...
                'DisplayName', 'trap, z_{meso}');
        else
            d3 = plot(y(mask), x(mask), 'o', 'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', mycolours(2,:), 'LineWidth', 0.5,...
                'DisplayName', 'trap, z_{meso}');
        end
        max3 = max(y(mask));
        hold on
        % zbathy - sediment traps
        mask = strcmp(g,'zbathy');
        if (sum(mask) == 0)
            d4 = plot(NaN, NaN, 'o', 'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', mycolours(1,:), 'LineWidth', 0.5,...
                'DisplayName', 'trap, z_{bathy}');
        else
            d4 = plot(y(mask), x(mask), 'o', 'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', mycolours(1,:), 'LineWidth', 0.5,...
                'DisplayName', 'trap, z_{bathy}');
        end
        max4 = max(y(mask));
        hold off

        box on
        
        if (iLoc == 1) % EqPac
            xlim([0 500])
            xTickValues = 0:200:400;
        elseif (iLoc == 2) % OSP
            xlim([0 200])
            xTickValues = 0:75:150;       
        elseif (iLoc == 3) % PAP-SO
            xlim([0 380])
            xTickValues = 0:150:300;
        elseif (iLoc == 4) % BATS/OFP
            xlim([0 350])
            xTickValues = 0:150:300;
        elseif (iLoc == 5) % HOT/ALOHA
            xlim([0 100])
            xTickValues = 0:50:100;
        elseif (iLoc == 6) % HAUSGARTEN
            xlim([0 70])
            xTickValues = 0:25:50;
        end
        xticks(xTickValues)
        xticklabels(xTickValues)

        ylim([15 4000])
        yticks([20 100 1000 4000])
        set(gca,'Yscale','log')
        axh = gca;
        axh.YAxis.TickDirection = 'out';
        axh.TickLength = [0.03, 0.03]; % make tick marks longer

        if (iMonth == 1 || iMonth == 4 || iMonth ==7 || iMonth == 10)
            yticklabels({'20','100','1000','4000'})
        else
            yticklabels([])
        end

        set(gca,'YDir','Reverse','XAxisLocation','Bottom','xlabel',[],'ylabel',[])

        title(monthLabel(iMonth),'FontSize',14)

    end % iMonth

    % Shift all plots a little bit to the right and up
    ax(1).pos(1) = ax(1).pos(1)+0.040; ax(1).pos(2) = ax(1).pos(2)+0.040; 
    ax(2).pos(1) = ax(2).pos(1)+0.025; ax(2).pos(2) = ax(2).pos(2)+0.040; 
    ax(3).pos(1) = ax(3).pos(1)+0.010; ax(3).pos(2) = ax(3).pos(2)+0.040;  
    ax(4).pos(1) = ax(4).pos(1)+0.040; ax(4).pos(2) = ax(4).pos(2)+0.040; 
    ax(5).pos(1) = ax(5).pos(1)+0.025; ax(5).pos(2) = ax(5).pos(2)+0.040; 
    ax(6).pos(1) = ax(6).pos(1)+0.010; ax(6).pos(2) = ax(6).pos(2)+0.040; 
    ax(7).pos(1) = ax(7).pos(1)+0.040; ax(7).pos(2) = ax(7).pos(2)+0.040; 
    ax(8).pos(1) = ax(8).pos(1)+0.025; ax(8).pos(2) = ax(8).pos(2)+0.040; 
    ax(9).pos(1) = ax(9).pos(1)+0.010; ax(9).pos(2) = ax(9).pos(2)+0.040; 
    ax(10).pos(1) = ax(10).pos(1)+0.040; ax(10).pos(2) = ax(10).pos(2)+0.040; 
    ax(11).pos(1) = ax(11).pos(1)+0.025; ax(11).pos(2) = ax(11).pos(2)+0.040; 
    ax(12).pos(1) = ax(12).pos(1)+0.010; ax(12).pos(2) = ax(12).pos(2)+0.040; 
    for iMonth = 1:12
        set(haxis(iMonth),'Position',ax(iMonth).pos) 
    end

    % Some touches to the legend
    lg = legend([d1 d2 d3 d4]);
    lg.Position(1) = 0.35; lg.Position(2) = -0.005;
    lg.Orientation = 'horizontal';
    lg.ItemTokenSize = [20,50];
    lg.FontSize = 11;
    lg.NumColumns = 2;
    lg.Box = 'on';

    % Give common xlabel, ylabel and title to your figure
    % Create a new axis
    a = axes;
    t = title(STATION_NAMES(iLoc),'FontSize',16);
    xl = xlabel('POC flux (mg C m^{-2} d^{-1})','FontSize',16);
    yl = ylabel('Depth (m)','FontSize',16);
    % Specify visibility of the current axis as 'off'
    a.Visible = 'off';
    % Specify visibility of Title, XLabel, and YLabel as 'on'
    t.Visible = 'on';
    xl.Visible = 'on';
    yl.Visible = 'on';
    yl.Position(1) = yl.Position(1) - 0.005; yl.Position(2) = yl.Position(2) + 0.02; 
    xl.Position(1) = 0.5; xl.Position(2) = xl.Position(2) + 0.065;
    t.Position(1) = t.Position(1); t.Position(2) = t.Position(2) + 0.035;

    saveFigure(strcat('compilation_pocflux_',STATION_TAGS{iLoc}))

end % iLoc

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - PLOT FIGURE S1 (NO. ENTRIES BY MONTH, LOCATION AND DEPTH
% HORIZON)
% -------------------------------------------------------------------------

theNumberOfDataPoints = classicMonthlyDhN(:,:,:); % nDepths x 12 months x nLocs
theNumberOfDataPoints_permutted = permute(theNumberOfDataPoints, [2 3 1]); % 12 months x nLocs x nDepths

% Swap locations: currently, the order is (1) EqPac, (2) OSP, (3) PAP-SO, 
% (4) BATS/OFP, (5) HOT/ALOHA and (6) HAUSGARTEN, and the desired order is 
% (1) HOT/ALOHA, (2) BATS/OFP, (3) EqPac, (4) PAP-SO, (5) OSP and (6) HAUSGARTEN.
theNumberOfDataPoints_swapped = theNumberOfDataPoints_permutted;
theNumberOfDataPoints_swapped(:, [3, 5, 4, 2, 1, 6], :) = theNumberOfDataPoints_permutted(:, [1, 2, 3, 4, 5, 6], :);

% .........................................................................

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.42 0.20],'Color','w')
axh = axes('Position', [0.10 0.27 0.82 0.62]);

y = squeeze(theNumberOfDataPoints_swapped(:,:,:));
yy = flipdim(y,3); % to have bathypelagic at the bottom of the plot instead of at the top
h = plotBarStackGroups(yy, monthLabel); % plot groups of stacked bars

% Change the colors of each bar segment
coloursBarSegments = parula(size(h,2));
coloursBarSegments = repelem(coloursBarSegments,size(h,1),1); 
coloursBarSegments = mat2cell(coloursBarSegments,ones(size(coloursBarSegments,1),1),3);
set(h,{'FaceColor'},coloursBarSegments)

ylim([0 200])
yl = ylabel('Number of data points');
yl.Position(1) = yl.Position(1) - 0.5;
box on
title('POC flux')

% Legend
lg = legend('Near seafloor','Base of mesopelagic','Base of euphotic');
lg.Position(1) = 0.40; lg.Position(2) = -0.023;
lg.Orientation = 'horizontal';
lg.FontSize = 12;
set(lg,'Box','off') 

set(gca, 'FontSize', 12); 

saveFigure('compilation_numberdatapoints')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - PLOT PERCENTAGE UNCERTAINTY BY MONTH, LOCATION AND DEPTH
% HORIZON
% -------------------------------------------------------------------------

fe = (classicMonthlyDhErrTot./classicMonthlyDhAvg).*100;
theFractionalError = fe(:,:,:); % nDepths x 12 months x nLocs
theFractionalError_permutted = permute(theFractionalError, [2 3 1]); % 12 months x nLocs x nDepths

% Swap locations: currently, the order is (1) EqPac, (2) OSP, (3) PAP-SO 
% and (4) BATS/OFP and the desired order is (1) BATS/OFP, (2) EqPac, 
% (3) PAP-SO and (4) OSP.
theFractionalError_swapped = theFractionalError_permutted;
theFractionalError_swapped(:, [3, 5, 4, 2, 1, 6], :) = theFractionalError_permutted(:, [1, 2, 3, 4, 5, 6], :);

% .........................................................................

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.42 0.20],'Color','w') 
axh = axes('Position', [0.10 0.27 0.82 0.62]);

y = squeeze(theFractionalError_swapped(:,:,:));
yy = flipdim(y,3); % to have bathypelagic at the bottom of the plot instead of at the top
h = plotBarStackGroups(yy, monthLabel); % plot groups of stacked bars

% Change the colors of each bar segment
coloursBarSegments = parula(size(h,2)); 
coloursBarSegments = repelem(coloursBarSegments,size(h,1),1); 
coloursBarSegments = mat2cell(coloursBarSegments,ones(size(coloursBarSegments,1),1),3);
set(h,{'FaceColor'},coloursBarSegments)

ylim([0 150])
yl = ylabel('% uncertainty');
yl.Position(1) = yl.Position(1) - 0.5;
box on
title('POC flux')

% Legend
lg = legend('Near seafloor','Base of mesopelagic','Base of euphotic');
lg.Position(1) = 0.40; lg.Position(2) = -0.023;
lg.Orientation = 'horizontal';
lg.FontSize = 12; 
set(lg,'Box','off') 

set(gca,'FontSize',12)

saveFigure('compilation_percentageuncertainty')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - PLOT FIGURE S2 (MONTHLY FLUXES AND THEIR ERROR BY LOCATION
% AND DEPTH HORIZON)
% -------------------------------------------------------------------------

err = squeeze(classicMonthlyDhErrTot(:,:,:)); % nDepths x 12 months x nLocs
vals = squeeze(classicMonthlyDhAvg(:,:,:)); % nDepths x 12 months x nLocs
parulaColours = flipud(parula(3));

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.70 0.50],'Color','w') 
haxis = zeros(3,NUM_LOCS);
iSubplot = 0;

for iDepthLayer = 1:3

    for iStation = 1:NUM_LOCS
        
        % Re-order
        switch iStation
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
        
        iSubplot = iSubplot + 1;

        haxis(iSubplot) = subaxis(3,NUM_LOCS,iSubplot,'Spacing',0.01,'Padding',0.01,'Margin', 0.07);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');

        if (iSubplot >= 1 && iSubplot <= 6)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.04;
        elseif (iSubplot >= 6 && iSubplot <= 12)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.01; 
         elseif (iSubplot >= 12 && iSubplot <= 18)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2) - 0.02;  
        end
        set(haxis(iSubplot),'Position',ax(iSubplot).pos)

        hbar = bar(haxis(iSubplot),(1:12),squeeze(vals(iDepthLayer,:,iLoc)),...
            'BarWidth',0.75,'FaceColor','flat');
        hbar.CData(:,:) = repmat(parulaColours(iDepthLayer,:),[12 1]);
        hold on
        
        her = errorbar(haxis(iSubplot),(1:12),squeeze(vals(iDepthLayer,:,iLoc)),...
            zeros(size(squeeze(vals(iDepthLayer,:,iLoc)))),squeeze(err(iDepthLayer,:,iLoc)));    
        her.Color = [0 0 0];                            
        her.LineStyle = 'none'; 
        hold off
        
        % Euphotic
        if (iSubplot >= 1 && iSubplot <= 6)
            ylim([0 22])
            yticks([0,5,10,15,20])
            yticklabels({'0','5','10','15','20'});
            ytickformat('%.0f')
        % Mesopelagic
        elseif (iSubplot > 6 && iSubplot <= 12)
            ylim([0 2.2])
            yticks([0,0.5,1,1.5,2])
            yticklabels({'0','0.5','1','1.5','2'});
            ytickformat('%.1f')
        % Bathypelagic
        else
            ylim([0 2.2])
            yticks([0,0.5,1,1.5,2])
            yticklabels({'0','0.5','1','1.5','2'});
            ytickformat('%.1f')
        end
        
        xlim([0.5 12+0.5])
        xticks(1:12);
%         if (iMonth == 11 || iMonth == 12)
%             xlabel('Year');
%             xticklabels(yearVectorData);
%             xtickangle(90); 
%         else
%             xticklabels([]);
%         end
        xticklabels(monthLabel);
        xtickangle(90); 
        axh = gca;
        axh.XAxis.FontSize = 9; 
        
        if (iSubplot >= 1 && iSubplot <= 6)
            tl = title(STATION_NAMES(iLoc),'FontSize',14);
            tl.Visible = 'on';
            tl.Position(2) = tl.Position(2) + 0.20;
        end
%         if (iSubplot == 4)
%             yl = ylabel('POC flux (mg C m^{-2} d^{-1})','FontSize',14);
%             yl.Visible = 'on';
%             yl.Position(1) = yl.Position(1) - 1;
%         end
        
        grid on;
        axh.XGrid = 'off';
        axh.YGrid = 'on';

    end
    
end

saveFigure('compilation_flux_by_month_and_station')

