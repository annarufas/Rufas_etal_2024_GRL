
% ======================================================================= %
%                                                                         %
% This script plots the POC flux data derived from the UVP5 dataset. It   %
% has 3 sections:                                                         %
%   Section 1 - Presets.                                                  %
%   Section 2 - Plot number of casts.                                     %
%   Section 3 - Plot Figure S4 (the UVP5-derived POC fluxes compared to   %
%               the trap and radionuclide-derived measurements).          %
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
addpath(genpath('./resources/external/'));
addpath(genpath('./resources/internal/'));

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

% Filename declarations 
filenameInputUvpProcessedDataset45sc = 'pocflux_bisson_45sc.mat';
filenameInputPocFluxCompilation      = 'pocflux_compilation.mat';
filenameInputTimeseriesInformation   = 'timeseries_station_information.mat';

% Load the UVP5 data
load(fullfile('.','data','processed','UVP5',filenameInputUvpProcessedDataset45sc),...
    'uvpFluxByCastAvg','uvpFluxByCastErr','uvpMonthlyFluxProfileAvg','ecotaxaDepths',...
    'castMonthlyDistrib')
maxNumCastsPerMonth = max(castMonthlyDistrib,[],'all');
nDepths = numel(ecotaxaDepths);

% Load the trap and radionuclide compilation
load(fullfile('.','data','processed',filenameInputPocFluxCompilation),...
    'classicRawProfileValues','classicRawProfileDataType','classicRawProfileDepths')

% Load station information
load(fullfile('.','data','processed',filenameInputTimeseriesInformation),...
    'STATION_NAMES','STATION_TAGS')
nLocs = length(STATION_NAMES);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - PLOT NUMBER OF CASTS
% -------------------------------------------------------------------------

monthsLabel = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
categoricalLabelMonths = categorical(monthsLabel);
categoricalLabelMonths = reordercats(categoricalLabelMonths,monthsLabel);

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.42 0.27],'Color','w') 
y = castMonthlyDistrib(:,:);
h = bar(categoricalLabelMonths,y);

% Change the colours of each bar
coloursBars = jet(size(h,2)); 
coloursBars = mat2cell(coloursBars,ones(size(coloursBars,1),1),3);
set(h,{'FaceColor'},coloursBars)

box on
yl = ylabel('Number of UVP5 casts','FontSize',12);
lg = legend(STATION_NAMES);
set(lg,'Box','off') 

saveFigure('uvp_numbercasts_45sizeclasses')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - PLOT FIGURE S4 (THE UVP5-DERIVED POC FLUXES COMPARED TO THE
% TRAP AND RADIONUCLIDE-DERIVED MEASUREMENTS)
% -------------------------------------------------------------------------

plotDepths = 47.5:5:2005;
idxPlotDepths = NaN(numel(plotDepths),1);
for iDepth = 1:numel(plotDepths)
    [~, idxPlotDepths(iDepth)] = min(abs(ecotaxaDepths - plotDepths(iDepth)));
end  

for iLoc = 1:nLocs
    
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.25 0.75],'Color','w')
    haxis = zeros(12,1);

    isCastPlotted = 0;
    
    for iMonth = 1:12

        haxis(iMonth) = subaxis(4,3,iMonth,'Spacing',0.055,'Padding',0,'Margin',0.13);
        ax(iMonth).pos = get(haxis(iMonth),'Position');
        
        % Calculate xlim
        maxAnnualUvp = max(uvpFluxByCastAvg(:,:,:,iLoc),[],'all','omitnan'); % mg C m-2 d-1
        maxAnnualClassic = max(classicRawProfileValues(:,:,iLoc),[],'all','omitnan'); % mg C m-2 d-1 

        maxMonthlyObs = max(classicRawProfileValues(:,iMonth,iLoc),[],'omitnan');
        if (isnan(maxMonthlyObs)) 
            maxMonthlyObs = 0;
        end
        maxMonthlyEst = max(uvpFluxByCastAvg(:,:,iMonth,iLoc),[],'all','omitnan');
        if (isnan(maxMonthlyEst)) 
            maxMonthlyEst = 0;
        end
        maxFlux = max(maxAnnualClassic,maxAnnualUvp);

        % Plot estimates from the UVP5: casts
        valse = squeeze(uvpFluxByCastAvg(:,idxPlotDepths,iMonth,iLoc));
        valse(valse==0) = NaN;
        valseresh = reshape(valse.', [], 1);
        ze = ecotaxaDepths(idxPlotDepths);
        zerep = repmat(ze,[maxNumCastsPerMonth 1]);
%         d1 = plot(haxis(iMonth),valseresh,zerep,'Color','r','LineWidth',0.4,...
%             'DisplayName','Estimated (individual casts)');
        % Replace red line by red scatters
        d1 = scatter(haxis(iMonth),valseresh,zerep,10,'r','filled',...
            'DisplayName','Estimated (individual casts)');
        
        hold on
                
        % Plot estimates from the UVP5: monthly mean
        valse = squeeze(uvpMonthlyFluxProfileAvg(idxPlotDepths,iMonth,iLoc,1));
        valse(valse==0) = NaN;
        ze = ecotaxaDepths(idxPlotDepths);
        d2 = plot(haxis(iMonth),valse,ze,'Color','k','LineWidth',1.5,...
            'DisplayName','Estimated (mean)');
        hold on
        
        % Plot measurements from sediment traps and radionuclides
        % (combined)
        
        % First - sediment trap
        mask = strcmp(squeeze(classicRawProfileDataType(:,iMonth,iLoc)),'trap');
        if (sum(mask) == 0)
            d3 = plot(NaN,NaN,'o','MarkerEdgeColor',[0.85 0.85 0.85],...
                'MarkerFaceColor',[0.85 0.85 0.85],'LineWidth',0.5,...
                'DisplayName','Observed (trap)');
        else
            d3 = plot(classicRawProfileValues(mask,iMonth,iLoc),...
                classicRawProfileDepths(mask,iMonth,iLoc),'o','MarkerEdgeColor',[0.85 0.85 0.85],...
                'MarkerFaceColor',[0.85 0.85 0.85],'LineWidth', 0.5,... 
                'DisplayName','Observed (trap)');
        end
        hold on
        
        % Second - radionuclides
        mask = strcmp(squeeze(classicRawProfileDataType(:,iMonth,iLoc)),'radionuclide');
        if (sum(mask) == 0)
            d4 = plot(NaN,NaN,'+','MarkerEdgeColor',[0.6 0.6 0.6],...
                'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',1.5,...
                'DisplayName', 'Observed (radionuclide)');
        else
            d4 = plot(classicRawProfileValues(mask,iMonth,iLoc),...
                classicRawProfileDepths(mask,iMonth,iLoc),'+','MarkerEdgeColor',[0.6 0.6 0.6],...
                'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',1.5,...
                'DisplayName','Observed (radionuclide)');
        end
        hold off
        
        box on
        
        if (iLoc == 1) % EqPac
            xlim([0 500])
            xTickValues = 0:200:400;
        elseif (iLoc == 2) % OSP
            xlim([0 200])
            xTickValues = 0:75:150;
        elseif (iLoc == 3 && iMonth ~= 5) % PAP-SO
            xlim([0 380])
            xTickValues = 0:150:300;
        elseif (iLoc == 3 && iMonth == 5) % PAP-SO
            xlim([0 800])
            xTickValues = 0:400:800;
        elseif (iLoc == 4 && iMonth ~= 12) % BATS/OFP
            xlim([0 350])
            xTickValues = 0:150:300;
        elseif (iLoc == 4 && iMonth == 12) % BATS/OFP
            xlim([0 1000])
            xTickValues = 0:400:800;
        elseif (iLoc == 5 && iMonth ~= 8) % HOT/ALOHA
            xlim([0 100])
            xTickValues = 0:50:100;
        elseif (iLoc == 5 && iMonth == 8) % HOT/ALOHA
            xlim([0 200])
            xTickValues = 0:100:200;
        elseif (iLoc == 6 && iMonth == 6) % HAUSGARTEN
            xlim([0 150])
            xTickValues = 0:70:140;
        elseif (iLoc == 6 && iMonth == 7) % HAUSGARTEN
            xlim([0 900])
            xTickValues = 0:400:800;
        elseif (iLoc == 6 && (iMonth ~= 6 || iMonth ~= 7)) % HAUSGARTEN
            xlim([0 70])
            xTickValues = 0:25:50;
        end
        xticks(xTickValues)
        xticklabels(xTickValues)

        ylim([0 1510])
        yticks([0 500 1000 1500])
        if (iMonth == 1 || iMonth == 4 || iMonth ==7 || iMonth == 10)
            yticklabels({'0','500','1000','1500'})
        else
            yticklabels([])
        end
        axh = gca;
        axh.YAxis.TickDirection = 'out';
        axh.TickLength = [0.03, 0.03]; % make tick marks longer
        
        title(monthsLabel(iMonth),'FontSize',14)

        set(gca,'YDir','Reverse','XAxisLocation','Bottom','xlabel',[],'ylabel',[])
        
        if (castMonthlyDistrib(iMonth,iLoc) > 0 && isCastPlotted == 0 && iMonth > 2)
            lg = legend([d1 d2 d3 d4],'NumColumns',2);
            lg.Position(1) = 0.12; lg.Position(2) = 0.02;
            lg.ItemTokenSize = [15,1];
            lg.FontSize = 11;
            set(lg,'Box','off') 
            isCastPlotted = 1;
        end

    end % iMonth
    
    for iMonth = 1:12
        if mod(iMonth, 3) == 1
            ax(iMonth).pos(1) = ax(iMonth).pos(1) + 0.040;
        elseif mod(iMonth, 3) == 2
            ax(iMonth).pos(1) = ax(iMonth).pos(1) + 0.025;
        else
            ax(iMonth).pos(1) = ax(iMonth).pos(1) + 0.010;
        end
        ax(iMonth).pos(2) = ax(iMonth).pos(2) + 0.040;
        set(haxis(iMonth),'Position',ax(iMonth).pos) 
    end 
    
    % Give common xlabel, ylabel and title to your figure
    % Create a new axis
    a = axes;
    t = title(STATION_NAMES(iLoc),'FontSize',16);
    xl = xlabel('POC flux (mg C m^{-2} d^{-1})','FontSize',16); % ,'FontWeight','bold'
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

    saveFigure(strcat('uvp_pocflux_45sizeclasses_',STATION_TAGS{iLoc}))
  
end % iLoc 
