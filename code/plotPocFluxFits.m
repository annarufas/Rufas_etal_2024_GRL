
% ======================================================================= %
%                                                                         %
% This script plots annual POC flux profiles from our sediment trap and   %
% radionuclide compilationas as well as the Martin's b fit for all extant %
% fit scenarios found in the literature (Figure S3).                      %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 20 Nov 2024                                   %
%                                                                         %
% ======================================================================= %

close all; clear all; clc
addpath(genpath('./data/interim/'));
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
filenameInputFluxCompilation       = 'pocflux_compilation.mat';
filenameInputTimeseriesInformation = 'timeseries_station_information.mat';
filenameInputMetricsArrayChoices   = 'fitmetrics_all_combinations.mat';

% Load the data
load(fullfile('.','data','processed',filenameInputFluxCompilation),...
    'classicAnnualProfileAvg','classicAnnualProfileErrTot',...
    'classicAnnualProfileN','classicAnnualProfileDepths')

% Load information on stations
load(fullfile('.','data','processed',filenameInputTimeseriesInformation),...
    'STATION_NAMES','qZeuMonthly','LOC_DEPTH_HORIZONS','MAX_ZEU')
nLocs = length(STATION_NAMES);

% Define possible values for each choice
isMeansOfMeansOptions   = [0, 1];
isLogTransformedOptions = [0, 1];
isFluxNormalisedOptions = [0, 1];
choiceZrefOptions       = [1, 2, 3]; % 1=closest value to 100, 2=zeu, 3=inflexion point

% Indexes to locations
iE = 1; % EqPac
iO = 2; % OSP
iP = 3; % PAP-SO
iB = 4; % BATS/OFP
iHo = 5; % HOT/ALOHA
iHa = 6; % HAUSGARTEN

% Canonical Martin's b
canonicalMartinb = 0.858;

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - PLOT ANNUAL MEAN DATA AND FITS
% -------------------------------------------------------------------------

load(fullfile('.','data','processed',filenameInputMetricsArrayChoices),...
    'martinbAnnual')

% Calculate annual euphotic layer depth
qZeuAnnual = mean(qZeuMonthly,1,'omitnan');
                                           
colourFeatures = brewermap(3,'*Set1');

% Loop through each combination of options  
for isLogTransformed = isLogTransformedOptions
    j = find(isLogTransformedOptions == isLogTransformed);

    for isFluxNormalised = isFluxNormalisedOptions
        k = find(isFluxNormalisedOptions == isFluxNormalised);

        for choiceZref = choiceZrefOptions

            figure()
            set(gcf,'Units','Normalized','Position',[0.01 0.05 0.58 0.28],'Color','w')
            haxis = zeros(nLocs,1);

            for iSubplot = 1:nLocs

                haxis(iSubplot) = subaxis(1,nLocs,iSubplot,'Spacing',0.03,'Padding',0.0,'Margin',0.11);
                ax(iSubplot).pos = get(haxis(iSubplot),'Position');

                % Shift all plots a little bit to the left and up
                ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.040;
                ax(iSubplot).pos(2) = ax(iSubplot).pos(2)+0.040;  
                set(haxis(iSubplot),'Position',ax(iSubplot).pos) 

                % Rearrange the order of locations: subtropics - tropics - subpolar
                if (iSubplot == 1)
                    iLoc = iHo; % HOT/ALOHA
                elseif (iSubplot == 2)
                    iLoc = iB;  % BATS/OFP
                elseif (iSubplot == 3)        
                    iLoc = iE;  % EqPac
                elseif (iSubplot == 4)         
                    iLoc = iP;  % PAP-SO
                elseif (iSubplot == 5)
                    iLoc = iO;  % OSP
                elseif (iSubplot == 6)
                    iLoc = iHa; % HAUSGARTEN
                end

                depthBounds = [max(LOC_DEPTH_HORIZONS(:,iLoc,1,1)), max(LOC_DEPTH_HORIZONS(:,iLoc,1,2)), max(LOC_DEPTH_HORIZONS(:,iLoc,1,3));
                               max(LOC_DEPTH_HORIZONS(:,iLoc,2,1)), max(LOC_DEPTH_HORIZONS(:,iLoc,2,2)), max(LOC_DEPTH_HORIZONS(:,iLoc,2,3))];

                % Extract data (as it is, not the 5-m bin-averaged)
                idxValidData = find(~isnan(classicAnnualProfileDepths(:,iLoc)));
                validDepths = classicAnnualProfileDepths(idxValidData,iLoc);
                validFluxes = classicAnnualProfileAvg(idxValidData,iLoc);
                validErrs = classicAnnualProfileErrTot(idxValidData,iLoc);
                [selectedDepths,selectedFluxes,selectedErrs] = extractDataFromZrefToEnd(...
                    validDepths,validFluxes,validErrs,choiceZref,depthBounds,MAX_ZEU);

                % Extract estimated fit for 'means of means' 
                estimatedMartinbMeansOfMeans = martinbAnnual(iLoc,1,j,k,choiceZref);

                % Extract estimated fit for 'annual fit' 
                estimatedMartinbAnnualFit = martinbAnnual(iLoc,2,j,k,choiceZref);

                % Plot annual mean of observations
                plot(haxis(iSubplot),validFluxes,validDepths,...
                    'o','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1.5);
                hold on;

%                 % Calculate the fit line for Martin's canonical fit and plot it
%                 [fluxFit,depthFit] = calculateFluxFit(selectedFluxes,selectedDepths,canonicalMartinb);
%                 plot(haxis(iSubplot),fluxFit,depthFit,'Color',colourFeatures(1,:),'LineWidth',2);
%                 hold on;

                % Calculate the fit line for estimated b for 'means of means' and plot it
                [fluxFit,depthFit] = calculateFluxFit(selectedFluxes,selectedDepths,estimatedMartinbMeansOfMeans);
                plot(haxis(iSubplot),fluxFit,depthFit,'Color',colourFeatures(2,:),'LineWidth',2);
                hold on;

                % Calculate the fit line for estimated b for 'annual fit' and plot it
                [fluxFit,depthFit] = calculateFluxFit(selectedFluxes,selectedDepths,estimatedMartinbAnnualFit);
                plot(haxis(iSubplot),fluxFit,depthFit,'Color',colourFeatures(3,:),'LineWidth',2);
                hold on;

                % Add reference depth line
                yline(selectedDepths(1),'--k','LineWidth',1);
                hold off;

                if (iLoc == 1) % EqPac
                    xlim([0 320])
                    xTickValues = 0:150:300;
                elseif (iLoc == 2) % OSP
                    xlim([0 220])
                    xTickValues = 0:100:200;       
                elseif (iLoc == 3) % PAP-SO
                    xlim([0 320])
                    xTickValues = 0:150:300;
                elseif (iLoc == 4) % BATS/OFP
                    xlim([0 55])
                    xTickValues = 0:25:50;
                elseif (iLoc == 5) % HOT/ALOHA
                    xlim([0 55])
                    xTickValues = 0:25:50;
                elseif (iLoc == 6) % HAUSGARTEN
                    xlim([0 25])
                    xTickValues = 0:10:20;
                end
                xticks(xTickValues)
                xticklabels(xTickValues)

                ylim([15 5000])
                yticks([20 100 1000 4000])
                set(gca,'Yscale','log')
                axh = gca;
                axh.YAxis.TickDirection = 'out';
                axh.TickLength = [0.03, 0.03]; % make tick marks longer

                yticklabels({'20','100','1000','4000'})
                set(gca,'YDir','Reverse','XAxisLocation','Bottom')

                if (iSubplot == nLocs)
                    lg = legend('Observed data',...
                                'Estimated b fit MM',...
                                'Estimated b fit AF',...
                                'Reference depth',...
                        'Location','eastoutside');
                    lg.Position(1) = 0.86; lg.Position(2) = 0.69;
                    lg.Orientation = 'vertical';
                    lg.FontSize = 12; 
                    lg.ItemTokenSize = [12,1];
                    set(lg,'Box','off')   
                end

                if (iSubplot == 1)
                    ylabel('Depth (m)','FontSize',12,'Units','normalized','Position',[-0.4, 0.5, 0])
                else
                    set(gca,'ylabel',[])
                    set(gca,'yticklabel',[])
                end

                if (iSubplot == 3)
                    xlabel('POC flux (mg C m^{-2} d^{-1})','FontSize',12,'Units','normalized','Position',[1.2, -0.10, 0])
                else
                    set(gca,'xlabel',[])
                end

                title(STATION_NAMES(iLoc),'FontSize',14)

            end % iSubplot

            % Create figure name using the 'filenameFitMetricsOutput'
            filenameFitMetricsOutput = constructFilenameFitMetrics(...
                1,isLogTransformed,isFluxNormalised,choiceZref);

            % Split the filename by underscores
            parts = split(filenameFitMetricsOutput, '_');
            part1 = parts{4}; % isLogTransformed
            part2 = parts{5}; % isFluxNormalised
            [part3,~] = strtok(parts{6},'.'); % choiceZref

            saveFigure(strcat('fitcomp_',part1,'_',part2,'_',part3))

        end % choiceZref
    end % isFluxNormalised
end % isLogTransformed

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS USED IN THIS SCRIPT
% -------------------------------------------------------------------------

function [fluxFit,depthRange] = calculateFluxFit(fluxValues,depthValues,b)

    % Parameters for the fit line
    z0 = depthValues(1); % initial depth (reference depth)
    F0 = fluxValues(1); % initial flux (at z0)

    % Generate the depth range of interest for the fit line
    depthRange = linspace(min(depthValues), max(depthValues), 100); 

    % Calculated flux values over the depth range
    fluxFit = F0.*(depthRange./z0).^(-b);

end
