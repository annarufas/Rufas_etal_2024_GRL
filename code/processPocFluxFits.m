
% ======================================================================= %
%                                                                         %
% This script calculates the Martin's b and remineralisation length scale %
% (z*) fits for all extant fit scenarios found in the literature for our  %
% sediment trap and radionuclide compilation. Additionally, the script    %
% writes the resulting fits to a .csv file (Dataset S2) and performs      %
% Principal Component Analysis (PCA) to explore potential groupings or    %
% patterns in the data.                                                   %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 6 Nov 2024                                    %
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

% Do we need to run the fit calculations?
areFitsComputed = 1;

% Filename declarations 
filenameInputTimeseriesInformation = 'timeseries_station_information.mat';
filenameInputFitMetricsArrayChoices = 'fitmetrics_all_combinations.mat';
filenameOutputCsvTable = 'dataset_s1_figureS4.csv';

% Load information on stations
load(fullfile('.','data','processed',filenameInputTimeseriesInformation),...
    'NUM_LOCS','STATION_NAMES','LOC_LATS','LOC_LONS')

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

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - EUPHOTIC LAYER DEPTH (NEEDED TO CALCULATE ZREF)
% -------------------------------------------------------------------------

% Load the global-ocean euphotic layer depth product calculated from CMEMS
% light attenuation coefficient (kd) product and extract data at our locations
    
if ~areFitsComputed
    
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
    qZeuAnnual = mean(qZeuMonthly,2,'omitnan');

end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - CALCULATE DIFFERENT FITS
% -------------------------------------------------------------------------

% Running this section might take up to 12 h

if ~areFitsComputed
    
    % Initialise output arrays
    martinbMonthly = NaN(NUM_LOCS,12,numel(isMeansOfMeansOptions),numel(isLogTransformedOptions),numel(isFluxNormalisedOptions),numel(choiceZrefOptions));
    zstarMonthly   = NaN(size(martinbMonthly)); 
    martinbAnnual  = NaN(NUM_LOCS,numel(isMeansOfMeansOptions),numel(isLogTransformedOptions),numel(isFluxNormalisedOptions),numel(choiceZrefOptions));
    zstarAnnual    = NaN(size(martinbAnnual));

    % Loop through each combination of options
    for isMeansOfMeans = isMeansOfMeansOptions
        i = find(isMeansOfMeansOptions == isMeansOfMeans);

        for isLogTransformed = isLogTransformedOptions
            j = find(isLogTransformedOptions == isLogTransformed);

            for isFluxNormalised = isFluxNormalisedOptions
                k = find(isFluxNormalisedOptions == isFluxNormalised);

                for choiceZref = choiceZrefOptions

                    % Display the current combination
                    fprintf('Running with isMeansOfMeans=%d, isLogTransformed=%d, isFluxNormalised=%d, choiceZref=%d\n', ...
                        isMeansOfMeans, isLogTransformed, isFluxNormalised, choiceZref);

                    % Call the function with the current combination
                    calculateBcpMetricsFromTrapAndRadionuclide(isMeansOfMeans, ...
                        isLogTransformed,isFluxNormalised,choiceZref);

                    % Load the data
                    filenameFitMetricsOutput = constructFilenameFitMetrics(...
                        isMeansOfMeans,isLogTransformed,isFluxNormalised,choiceZref);
                    data = load(filenameFitMetricsOutput,...
                        'martinbMonthly','zstarMonthly','martinbAnnual','zstarAnnual');

                    % Offload into output arrays
                    if isMeansOfMeans  
                        martinbMonthly(:,:,i,j,k,choiceZref) = data.martinbMonthly(:,:,1);
                        zstarMonthly(:,:,i,j,k,choiceZref) = data.zstarMonthly(:,:,1);
                    end
                    martinbAnnual(:,i,j,k,choiceZref) = data.martinbAnnual(:,1);
                    zstarAnnual(:,i,j,k,choiceZref) = data.zstarAnnual(:,1);

                end % choiceZref
            end % isFluxNormalised
        end % isLogTransformed
    end % isMeansOfMeans

    save(fullfile('.','data','processed',filenameInputFitMetricsArrayChoices),...
        'martinbMonthly','zstarMonthly','martinbAnnual','zstarAnnual','qZeuAnnual')

end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - WRITE OUT THE FITS DATA 5D ARRAY INTO A .CSV FILE (DATASET S1)
% -------------------------------------------------------------------------

load(fullfile('.','data','processed',filenameInputFitMetricsArrayChoices),...
    'martinbAnnual','zstarAnnual')

% I need to flatten or reshape the array into a 2D matrix 

% Initialise the output array
outputArray = NaN(NUM_LOCS,(2*2*2*2*3));

%                                                   Annual Martin's b fit                                                   Annual zstar fit
%           ------------------------------------------------------------------------------------------------------- ----------------------------------
%                                    Mean of monthly fits                           Fit to annual POC flux profiles        
%           ----------------------------------------------------------------------- ------------------------------- ----------------- ----------------
%                        not log-log                           log-log                         ...    ...               ...    ...       ...    ...
%           ----------------------------------- ----------------------------------- ------------------------------- ----------------------------------
%           not norm.  norm.  not norm.  norm.  not norm.  norm.  not norm.  norm.        ... ... ... ... ...        ... ... ... ...  ... ... ... ...
%           -------- -------- -------- -------- -------- -------- -------- -------- ------------------------------- ----------------------------------
%           z1 z2 z3 z1 z2 z3 z1 z2 z3 z1 z2 z3 z1 z2 z3 z1 z2 z3 z1 z2 z3 z1 z2 z3       . . . . . . . . . . .      ... ... ... ...  ... ... ... ...
%           -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- ------------------------------- ----------------------------------
% HOT/ALOHA
% BATS/OFP
% EqPac
% ...
% ...

funb = @(x) sprintf('%.2f', x); % decimal digits for Martin's b
funz = @(x) sprintf('%3.0f', x); % decimal digits for z*

iRow = 0; % index for the output array
iColumn = 0; % index for the output array

for iLoc = [iHo,iB,iE,iP,iO,iHa]
    iRow = iRow + 1;
    iColumn = 0;
    
    for iMetric = 1:2   
        for isMeansOfMeans = isMeansOfMeansOptions
            i = find(isMeansOfMeansOptions == isMeansOfMeans);

            for isLogTransformed = isLogTransformedOptions
                j = find(isLogTransformedOptions == isLogTransformed);

                for isFluxNormalised = isFluxNormalisedOptions
                    k = find(isFluxNormalisedOptions == isFluxNormalised);

                    for choiceZref = choiceZrefOptions
                        iColumn = iColumn + 1;
                        
                        if (iMetric == 1) % Martin's b
                            extract = martinbAnnual(iLoc,i,j,k,choiceZref);
                        elseif (iMetric == 2) % zstar
                            extract = zstarAnnual(iLoc,i,j,k,choiceZref);
                        end
                        
                        % Cap the number of decimal digits accordingly
                        if ~isnan(extract)
                            C = num2cell(extract);
                            if (iMetric == 1) % Martin's b
                                F = cellfun(funb,C,'UniformOutput',0);
                            elseif (iMetric == 2)
                                F = cellfun(funz,C,'UniformOutput',0);
                            end
                            extract = str2double(F);
                        end
                        outputArray(iRow,iColumn) = extract;

                    end % choiceZref
                end % isFluxNormalised
            end % isLogTransformed
        end % isMeansOfMeans 
    end % iMetric
end % iLoc

% Concatenate row names and numeric array horizontally...
horzConcatenatedArray = [{'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'}',num2cell(outputArray)];

% ... and now vertically with the headers
headerRowFirst = [{''},repmat({'Martin_b'},1,(2*2*2*3)),repmat({'z_star'},1,(2*2*2*3))];

headerRowSecond = [{''},repmat([repmat({'Mean of monthly fits'},1,(2*2*3)),repmat({'Fit to annual POC flux profiles'},1,(2*2*3))],1,2)];

headerRowThird = [{''},repmat([repmat({'not log-log'},1,(2*3)),repmat({'log-log'},1,(2*3))],1,4)];

headerRowFourth = [{''},repmat([repmat({'not normalised'},1,3),repmat({'normalised'},1,3)],1,8)];

headerRowFifth = ['Location',repmat({'100 m','zeu','inflexion'},1,(2*2*2*2))];

outputTable = cell2table([headerRowFirst;...
                          headerRowSecond;...
                          headerRowThird;...
                          headerRowFourth;...
                          headerRowFifth;...
                          horzConcatenatedArray]);

% Write table to CSV without column names
writetable(outputTable,fullfile('.','data','processed',filenameOutputCsvTable),...
    'WriteVariableNames',false,'Delimiter',',');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - PRINCIPAL COMPONENT ANALYSIS
% -------------------------------------------------------------------------

% Select Martin's b only
martinbData = outputArray(:,1:24);
martinbData = martinbData'; % transpose
martinbData(isnan(martinbData)) = 0; % PCA does not like NaN
martinbData = martinbData(:,1:5); % no HAUSGARTEN

% Create the labels vector with the different categories
groupLabels = {'MM, not log-log, not norm, 100 m',...
               'MM, not log-log, not norm, zeu',...
               'MM, not log-log, not norm, inflex',...
               'MM, not log-log, norm, 100 m',...
               'MM, not log-log, norm, zeu',...
               'MM, not log-log, norm, inflex',...
               'MM, log-log, not norm, 100 m',...
               'MM, log-log, not norm, zeu',...
               'MM, log-log, not norm, inflex',...
               'MM, log-log, norm, 100 m',...
               'MM, log-log, norm, zeu',...
               'MM, log-log, norm, inflex',...
               'AF, not log-log, not norm, 100 m',...
               'AF, not log-log, not norm, zeu',...
               'AF, not log-log, not norm, inflex',...
               'AF, not log-log, norm, 100 m',...
               'AF, not log-log, norm, zeu',...
               'AF, not log-log, norm, inflex',...
               'AF, log-log, not norm, 100 m',...
               'AF, log-log, not norm, zeu',...
               'AF, log-log, not norm, inflex',...
               'AF, log-log, norm, 100 m',...
               'AF, log-log, norm, zeu',...
               'AF, log-log, norm, inflex'}';

% Define the colors based on the categories
myColourMap = brewermap(8,'Spectral');

% Initialise the color and symbol arrays
colourArray = zeros(length(groupLabels),3);
symbolArray = cell(length(groupLabels),1); 

% Loop through the group labels to assign colors and symbols
for i = 1:length(groupLabels)
    label = groupLabels{i};
    
    % Determine the symbol based on the suffix (100 m, zeu, inflex)
    if contains(label, '100 m')
        symbolArray{i} = 'o'; 
    elseif contains(label, 'zeu')
        symbolArray{i} = '+';  
    elseif contains(label, 'inflex')
        symbolArray{i} = 'x';
    end
    
    % Determine the color based on the prefix (MM, AF) and middle part
    if contains(label, 'MM')
        if contains(label, 'not log-log, not norm')
            colourArray(i, :) = myColourMap(1,:);
        elseif contains(label, 'not log-log, norm')
            colourArray(i, :) = myColourMap(2,:); 
        elseif contains(label, 'log-log, not norm')
            colourArray(i, :) = myColourMap(3,:);
        elseif contains(label, 'log-log, norm')
            colourArray(i, :) = myColourMap(4,:); 
        end
    elseif contains(label, 'AF')
        if contains(label, 'not log-log, not norm')
            colourArray(i, :) = myColourMap(5,:);
        elseif contains(label, 'not log-log, norm')
            colourArray(i, :) = myColourMap(6,:); 
        elseif contains(label, 'log-log, not norm')
            colourArray(i, :) = myColourMap(7,:);
        elseif contains(label, 'log-log, norm')
            colourArray(i, :) = myColourMap(8,:);
        end
    end
end

% Perform PCA
[coeff,score,latent,~,explained] = pca(zscore(martinbData)); % standardise the data (z-score normalisation) before applying PCA

% Create a scatter plot of PC1 vs PC2
figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.40],'Color','w')
h = gscatter(score(:,1), score(:,2), groupLabels, colourArray, char(symbolArray));
for i = 1:length(h)
    h(i).MarkerSize = 12; 
    h(i).LineWidth = 2;    
end
xlabel(['PC1 (' num2str(explained(1), '%.0f') '% variance explained)'],'FontSize',12);
ylabel(['PC2 (' num2str(explained(2), '%.0f') '% variance explained)'],'FontSize',12);
title('Principal Component Analysis','FontSize',12);
lg = legend(groupLabels,'Location','eastoutside');
lg.Orientation = 'vertical';
lg.FontSize = 9; 
lg.Title.String = 'Categories';  
set(lg,'Box','off')   
grid on
exportgraphics(gcf,fullfile('.','figures','fitcomp_PCA.png'),'Resolution',600)

