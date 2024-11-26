
% ======================================================================= %
%                                                                         %
% This script reads in the POC flux compilation created for this study    %
% (Dataset S0), calculates monthly and annual averages and propagates     %
% error accordingly. The script has 9 sections:                           %
%   Section 1 - Presets.                                                  %
%   Section 2 - Load the dataset and manipulate the data array.           %
%   Section 3 - Get euphotic layer depth.                                 %
%   Section 4 - Bin data monthly by depth horizon and propagate error.    %
%   Section 5 - Bin data monthly by unique depth and propagate error.     %
%   Section 6 - Bin data annually by depth horizon and propagate error.   %
%   Section 7 - Bin data annually by unique depth and propagate error.    %
%   Section 8 - Calculate the number of data points based on various      % 
%               criteria.                                                 %
%   Section 9 - Save the data.                                            %              
%                                                                         %
%   This script uses these external functions:                            % 
%       cleverTimeInterpolation.m - custom function                       %
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
filenameInputFluxCompilation        = 'dataset_s0_trap_and_radionuclide_compilation.xlsx';
filenameOutputFluxCompilation       = 'pocflux_compilation.mat';
filenameOutputTimeseriesInformation = 'timeseries_station_information.mat';

% Define parameters
POCFLUX_RAND_ERR_FRAC = 0.30; % 30% (Buesseler et al. 2000, Buesseler et al. 2007, Stanley et al. 2004)
POCFLUX_SYS_ERR_FRAC = 0.10; % 10%, based on a literature review
OC_ERR_FRAC = 0.10; % 10%, based on McKinna et al. 2019
MAX_ZEU = 200; % m
MOLAR_MASS_CARBON = 12.011; % g mol-1
MAX_NUM_VALUES_PER_MONTH = 1000;
MAX_NUM_DEPTHS_PER_PROFILE = 100;

% Load the global-ocean euphotic layer depth product calculated from 
% CMEMS' kd used to calculate zeu
load(fullfile('.','data','interim','zeu_calculated_kdcmems_mldcmems_pointonepercentpar0.mat'),...
    'zeu','zeu_lat','zeu_lon')

% Station information
STATION_NAMES = {'EqPac','OSP','PAP-SO','BATS/OFP','HOT/ALOHA','HAUSGARTEN'}; 
STATION_TAGS = {'eqpac','osp','papso','batsofp','hotaloha','hausgarten'}; 
nLocs = length(STATION_NAMES);

% The raw POC flux data array should contain data for the following ocean
% sites, which have the following latitude and longitudes:

          % EqPac  OSP PAP-SO   BATS  HOT HAUSGARTEN
LOC_LATS = [   0,   50,    49,  31.6, 22.5,  79]; % HOT/ALOHA lat slightly modified to extract data from WOA 
LOC_LONS = [-140, -145, -16.5, -64.2, -158, 4.3];

% The following depth horizons defined on the data will be used for data 
% summaries
LOC_DEPTH_HORIZONS = zeros(12,nLocs,2,3); % 2 boundaries for each of the 3 depth horizons

% EqPac                              % OSP                                % PAP-SO                          
% Targets data at zeu                Targets data at zeu                  Targets data at zeu
LOC_DEPTH_HORIZONS(:,1,1,1) = NaN;   LOC_DEPTH_HORIZONS(:,2,1,1) = NaN;   LOC_DEPTH_HORIZONS(:,3,1,1) = NaN;   
LOC_DEPTH_HORIZONS(:,1,2,1) = NaN;   LOC_DEPTH_HORIZONS(:,2,2,1) = NaN;   LOC_DEPTH_HORIZONS(:,3,2,1) = NaN;   
% Targets the trap at 880 m          Targets the traps at 1000 & 1009 m   Targets the trap at 1000 m     
LOC_DEPTH_HORIZONS(:,1,1,2) = 870;   LOC_DEPTH_HORIZONS(:,2,1,2) = 990;   LOC_DEPTH_HORIZONS(:,3,1,2) = 990;  
LOC_DEPTH_HORIZONS(:,1,2,2) = 890;   LOC_DEPTH_HORIZONS(:,2,2,2) = 1010;  LOC_DEPTH_HORIZONS(:,3,2,2) = 1010; 
% Targets the trap at 3618 m         Targets the traps at 3800 & 3805 m   Targets the trap at 4700 m
LOC_DEPTH_HORIZONS(:,1,1,3) = 3610;  LOC_DEPTH_HORIZONS(:,2,1,3) = 3790;  LOC_DEPTH_HORIZONS(:,3,1,3) = 4690; 
LOC_DEPTH_HORIZONS(:,1,2,3) = 3630;  LOC_DEPTH_HORIZONS(:,2,2,3) = 3810;  LOC_DEPTH_HORIZONS(:,3,2,3) = 4710; 

% BATS/OFP                           % HOT/ALOHA                          % HAUSGARTEN
% Targets data at zeu                Targets data at zeu                  Targets data at zeu 
LOC_DEPTH_HORIZONS(:,4,1,1) = NaN;   LOC_DEPTH_HORIZONS(:,5,1,1) = NaN;   LOC_DEPTH_HORIZONS(:,6,1,1) = NaN;
LOC_DEPTH_HORIZONS(:,4,2,1) = NaN;   LOC_DEPTH_HORIZONS(:,5,2,1) = NaN;   LOC_DEPTH_HORIZONS(:,6,2,1) = NaN;
% Targets the trap at 1500 m         Targets the trap at 1500 m           Targets the traps at 1225 & 1250 m
LOC_DEPTH_HORIZONS(:,4,1,2) = 1490;  LOC_DEPTH_HORIZONS(:,5,1,2) = 1490;  LOC_DEPTH_HORIZONS(:,6,1,2) = 1220;
LOC_DEPTH_HORIZONS(:,4,2,2) = 1510;  LOC_DEPTH_HORIZONS(:,5,2,2) = 1510;  LOC_DEPTH_HORIZONS(:,6,2,2) = 1260;
% Targets the trap at 3200 m         Targets the trap at 4000 m           Targets the traps at 2560 & 2618 m
LOC_DEPTH_HORIZONS(:,4,1,3) = 3190;  LOC_DEPTH_HORIZONS(:,5,1,3) = 3990;  LOC_DEPTH_HORIZONS(:,6,1,3) = 2550;
LOC_DEPTH_HORIZONS(:,4,2,3) = 3210;  LOC_DEPTH_HORIZONS(:,5,2,3) = 4010;  LOC_DEPTH_HORIZONS(:,6,2,3) = 2620;

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - LOAD THE DATASET AND MANIPULATE THE DATA ARRAY
% -------------------------------------------------------------------------

% Load the excel spreadhseet with flux data
opts = detectImportOptions(filenameInputFluxCompilation);
opts = setvartype(opts,{'POC_mmol_m2_d','POC_mg_m2_d','randerr_POC_mmol_m2_d','randerr_POC_mg_m2_d'},'double');
opts = setvartype(opts,{'deploymentDate','midDate','recoveryDate'},'datetime');
opts = setvartype(opts,{'method','source','originalUnits','comments'},'string');
D = readtable(filenameInputFluxCompilation,opts);

% Add a 'month' and a 'year' column
D.month = month(D.midDate);
D.year = year(D.midDate);

% Transform into categorical the variable 'month' and 'tag'
D.tag = categorical(D.tag);
D.month = categorical(D.month);

% Add a column to indicate the depth horizon of the data 
D.depthHorizon = cell(height(D),1);

% Make sure the timetable is sorted by 'tag' - ALPHABETICAL ORDER
D = sortrows(D,'tag');

% Add the sys error column
D.syserr_POC_mmol_m2_d = POCFLUX_SYS_ERR_FRAC .* D.POC_mmol_m2_d;
D.syserr_POC_mg_m2_d = MOLAR_MASS_CARBON .* D.syserr_POC_mmol_m2_d;

% Update the random error
D.randerr_POC_mmol_m2_d(isnan(D.randerr_POC_mmol_m2_d)) = POCFLUX_RAND_ERR_FRAC .* D.POC_mmol_m2_d(isnan(D.randerr_POC_mmol_m2_d));
D.randerr_POC_mg_m2_d = MOLAR_MASS_CARBON .* D.randerr_POC_mmol_m2_d;

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - GET EUPHOTIC LAYER DEPTH
% -------------------------------------------------------------------------

% Query points for interpolation
qLats = LOC_LATS;
qLons = LOC_LONS;

% Original data grid
[X,Y,T] = ndgrid(zeu_lat,zeu_lon,(1:12)');

% Interpolant 
F = griddedInterpolant(X, Y, T, zeu, 'linear'); 

% Extract data for the study locations (defined by qLats and qLons)
qZeuMonthly = NaN(12,nLocs);
for iLoc = 1:nLocs
    [qX,qY,qT] = ndgrid(qLats(iLoc),qLons(iLoc),(1:12)');
    qZeuMonthly_temp = F(qX,qY,qT);

    % Replace NaNs using an interpolation method
    if any(isnan(qZeuMonthly_temp))
        qZeuMonthly_temp(isnan(qZeuMonthly_temp)) = 0;
        qZeuMonthly(:,iLoc) = cleverTimeInterpolation(qZeuMonthly_temp,(1:12));
    else
        qZeuMonthly(:,iLoc) = qZeuMonthly_temp;
    end
    
    % Add information to LOC_DEPTH_HORIZONS array
    LOC_DEPTH_HORIZONS(:,iLoc,1,1) = qZeuMonthly(:,iLoc)... 
        - OC_ERR_FRAC.*qZeuMonthly(:,iLoc); % -10%
    LOC_DEPTH_HORIZONS(:,iLoc,2,1) = qZeuMonthly(:,iLoc)... 
        + OC_ERR_FRAC.*qZeuMonthly(:,iLoc); % +10%
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - BIN DATA MONTHLY BY DEPTH HORIZON AND PROPAGATE ERROR (THESE
% DATA WILL BE USED FOR SUMMARIES AND PLOTTING)
% -------------------------------------------------------------------------

% The variables that we want to extract
tagPocValues = 'POC_mmol_m2_d';
tagSysError = 'syserr_POC_mmol_m2_d';
tagRandError = 'randerr_POC_mmol_m2_d';
tagDepthHorizons = {'zeu','zmeso','zbathy'};

% Define output arrays
classicRawProfileValues_cell   = cell(12,nLocs); 
classicRawProfileErrRand_cell  = cell(12,nLocs); 
classicRawProfileErrSys_cell   = cell(12,nLocs);
classicRawProfileDepths_cell   = cell(12,nLocs); 
classicRawProfileDataType_cell = cell(12,nLocs); 

classicRawDhValues_cell   = cell(3,12,nLocs); 
classicRawDhDepths_cell   = cell(3,12,nLocs); 
classicRawDhDataType_cell = cell(3,12,nLocs); % sediment trap vs radionuclides
classicRawDhTag_cell      = cell(3,12,nLocs); % zeu, zmeso or zbathy

classicMonthlyDhAvg    = NaN(3,12,nLocs);   % mean
classicMonthlyDhN      = zeros(3,12,nLocs); % number of values
classicMonthlyDhErrTot = NaN(3,12,nLocs);   % total error

for iLoc = 1:nLocs
    currStatData = D(D.tag == STATION_NAMES{iLoc},:);

    % Locate the variables that are relevant
    colIdsAll = zeros(6,1);
    colIdsAll(1) = find(strcmpi(currStatData.Properties.VariableNames,tagPocValues));
    colIdsAll(2) = find(strcmpi(currStatData.Properties.VariableNames,tagSysError));
    colIdsAll(3) = find(strcmpi(currStatData.Properties.VariableNames,tagRandError));
    colIdsAll(4) = find(strcmpi(currStatData.Properties.VariableNames,'month'));
    colIdsAll(5) = find(strcmpi(currStatData.Properties.VariableNames,'depth'));
    colIdsAll(6) = find(strcmpi(currStatData.Properties.VariableNames,'method'));

    % Extract these variables
    currStatData = currStatData(:,colIdsAll');
            
    % Identify rows with any NaN values
    rowsWithNaN = any(ismissing(currStatData), 2);

    % Remove rows with NaN values
    currStatDataValid = currStatData(~rowsWithNaN, :);
    
    % (1) Pull all the data at that station and save it to output arrays
    for iMonth = 1:12
        currMonthData = currStatDataValid(currStatDataValid.month == num2str(iMonth),:);
        if (~isempty(currMonthData))
            classicRawProfileValues_cell{iMonth,iLoc} = num2str(currMonthData{:,1}'); 
            classicRawProfileErrSys_cell{iMonth,iLoc} = num2str(currMonthData{:,2}');
            classicRawProfileErrRand_cell{iMonth,iLoc} = num2str(currMonthData{:,3}'); 
            classicRawProfileDepths_cell{iMonth,iLoc} = num2str(currMonthData{:,5}');
            classicRawProfileDataType_cell{iMonth,iLoc} = strjoin(currMonthData{:,6}'); 
        end  
    end

    % (2) Pull data by month
    for iMonth = 1:12
        monthCondition = currStatDataValid.month == num2str(iMonth);
        currStatMonthData = currStatDataValid(monthCondition,:);
            
        % (3) Pull data by the depth horizon of interest
        for iDh = 1:3 
            depthCondition = currStatMonthData.depth > LOC_DEPTH_HORIZONS(iMonth,iLoc,1,iDh) & ...
                             currStatMonthData.depth < LOC_DEPTH_HORIZONS(iMonth,iLoc,2,iDh);
            
            if sum(depthCondition) ~= 0
                
                currStatMonthDhData = currStatMonthData(depthCondition,:);
                
            elseif sum(depthCondition) == 0
                
                % For iDh = 1, if there are no data between the depth layer 
                % boundaries defined in LOC_DEPTH_HORIZONS, find data at
                % the closest depth that is < 200 m
                if iDh == 1 && any(currStatMonthData.depth < MAX_ZEU)
                    absDiff = abs(currStatMonthData.depth - qZeuMonthly(iMonth,iLoc)); % calculate the absolute difference between each depth and locZeu
                    minDiff = min(absDiff); % find the minimum absolute difference
                    currStatMonthDhData = currStatMonthData(absDiff == minDiff,:); % extract data that match the minimum absolute difference
                
                    % Redefine LOC_DEPTH_HORIZONS based on the depth above
                    depthMatched = unique(currStatMonthDhData.depth);
                    LOC_DEPTH_HORIZONS(iMonth,iLoc,1,iDh) = depthMatched... 
                        - OC_ERR_FRAC.*depthMatched; % -10%
                    LOC_DEPTH_HORIZONS(iMonth,iLoc,2,iDh) = depthMatched...
                        + OC_ERR_FRAC.*depthMatched; % -10%
                
                else
                    currStatMonthDhData = [];
                end
            end
            
            % Proceed if there are data in the current month
            if ~isempty(currStatMonthDhData)

                % Add the depth horizon tag to currStatMonthDhData
                currStatMonthDhData.depthHorizon(:) = tagDepthHorizons(iDh);

                % Add the depth horizon tag to D array for later use
                idxRows = find(D.tag == STATION_NAMES{iLoc} &...
                        D.depth > LOC_DEPTH_HORIZONS(iMonth,iLoc,1,iDh) &...
                        D.depth < LOC_DEPTH_HORIZONS(iMonth,iLoc,2,iDh));
                D.depthHorizon(idxRows) = tagDepthHorizons(iDh);

                % Store values, their depth and their type
                classicRawDhValues_cell{iDh,iMonth,iLoc} = num2str(currStatMonthDhData{:,1}'); 
                classicRawDhDepths_cell{iDh,iMonth,iLoc} = num2str(currStatMonthDhData{:,5}');
                classicRawDhDataType_cell{iDh,iMonth,iLoc} = strjoin(currStatMonthDhData{:,6}'); 
                classicRawDhTag_cell{iDh,iMonth,iLoc} = strjoin(currStatMonthDhData{:,7}');

                % (4) Calculate the average & get the associated no. data points (N)
                classicMonthlyDhAvg(iDh,iMonth,iLoc) = mean(currStatMonthDhData{:,1}).*MOLAR_MASS_CARBON; % mmol C m-2 d-1 --> mg C m-2 d-1
                classicMonthlyDhN(iDh,iMonth,iLoc) = sum(currStatMonthDhData{:,1} >= 0,'omitnan');
                
                % (5) Calculate net error (propagate type A and type B
                % errors separately)
                vals = squeeze(currStatMonthDhData{:,1});
                errSys = squeeze(currStatMonthDhData{:,2});
                errRand = squeeze(currStatMonthDhData{:,3});
                errSys(isnan(errSys)) = 0;
                errRand(isnan(errRand)) = 0;

                % Error propagation of random errors (type A)
                [~,~,~,f_MID,f_UB,~,~] = worstcase(@(x) mean(x),vals,errRand);
                typeAuncertainty = f_UB - f_MID;

                % Error propagation of the systematic errors (type B)
                [~,~,~,f_MID,f_UB,~,~] = worstcase(@(x) mean(x),vals,errSys);
                typeBuncertainty = f_UB - f_MID;

                % Net error: sum of type A and type B errors in the quadrature
                classicMonthlyDhErrTot(iDh,iMonth,iLoc) =...
                    sqrt(typeAuncertainty^2 + typeBuncertainty^2).*MOLAR_MASS_CARBON; % mmol C m-2 d-1 --> mg C m-2 d-1

            end % are there data in the current month?
        end % iMonth    
    end % iDh
end % iLoc

% % Checks
% a_vals = squeeze(classicMonthlyDhAvg(1,:,:));
% a_err = squeeze(classicMonthlyDhErrTot(1,:,:));
% a_stdperc = (a_err./a_vals).*100;
% a_n = squeeze(classicMonthlyDhN(NUM_LOCS,:,:));

% Tide up the variable 'depthHorizon'. Find the rows that have not been 
% assigned into a depth horizon and assign them a 'NaN' string
isEmptyDhTag = cellfun(@isempty,D.depthHorizon);
iEmptyRows = find(isEmptyDhTag);
nEmptyRows = length(iEmptyRows);
D.depthHorizon(iEmptyRows) = repmat({'NaN'},nEmptyRows,1);
D.depthHorizon = categorical(D.depthHorizon);
TRAPRAD_TABLE = D;
DP = D(D.depthHorizon ~= 'NaN',:); % data processed

clear vals errRand errSys

% Save for use in other scripts
save(fullfile('.','data','processed',filenameOutputTimeseriesInformation),...
    'LOC_LATS','LOC_LONS','LOC_DEPTH_HORIZONS','STATION_NAMES','STATION_TAGS',...
    'MAX_NUM_VALUES_PER_MONTH','MAX_NUM_DEPTHS_PER_PROFILE','qZeuMonthly','MAX_ZEU')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - BIN DATA MONTHLY BY UNIQUE DEPTH AND PROPAGATE ERROR
% -------------------------------------------------------------------------

% Unfold the data stored in three cell arrays and store them into numerical 
% arrays.

classicRawProfileValues   = zeros(MAX_NUM_VALUES_PER_MONTH,12,nLocs);
classicRawProfileErrRand  = zeros(MAX_NUM_VALUES_PER_MONTH,12,nLocs);
classicRawProfileErrSys   = zeros(MAX_NUM_VALUES_PER_MONTH,12,nLocs);
classicRawProfileDepths   = zeros(MAX_NUM_VALUES_PER_MONTH,12,nLocs);
classicRawProfileDataType = cell(MAX_NUM_VALUES_PER_MONTH,12,nLocs);

for iLoc = 1:nLocs
    for iMonth = 1:12

        currMonthData     = zeros(MAX_NUM_VALUES_PER_MONTH,1);
        currMonthErrRand  = zeros(MAX_NUM_VALUES_PER_MONTH,1);
        currMonthErrSys   = zeros(MAX_NUM_VALUES_PER_MONTH,1);
        currMonthDepths   = zeros(MAX_NUM_VALUES_PER_MONTH,1);
        currMonthDataType = cell(MAX_NUM_VALUES_PER_MONTH,1);

        allmyvals = classicRawProfileValues_cell{iMonth,iLoc};
        allmyranderrs = classicRawProfileErrRand_cell{iMonth,iLoc};
        allmysyserrs = classicRawProfileErrSys_cell{iMonth,iLoc};
        allmydepths = classicRawProfileDepths_cell{iMonth,iLoc};
        allmydatatypes = classicRawProfileDataType_cell{iMonth,iLoc};
        
        if (~isempty(allmyvals))
            nDataPoints = length(str2num(allmyvals));
            currMonthData(1:nDataPoints) = str2num(allmyvals);
            currMonthErrSys(1:nDataPoints) = str2num(allmysyserrs);
            currMonthErrRand(1:nDataPoints) = str2num(allmyranderrs);
            thetypesall = textscan(allmydatatypes,'%s');
            for iDataPoint = 1:nDataPoints
                currMonthDataType{iDataPoint} = thetypesall{1}{iDataPoint};
            end
            currMonthDepths(1:nDataPoints) = str2num(allmydepths);
            
            % Correct for NaN
            theranderrs = currMonthErrRand(1:nDataPoints);
            theranderrs(isnan(theranderrs)) = 0;
            currMonthErrRand(1:nDataPoints) = theranderrs;
        else
            nDataPoints = 0;
        end
        
        [rowIdxs,~,yall] = find(currMonthData);
        xall = zeros(nDataPoints,1);
        gall = cell(nDataPoints,1);
        for iDataPoint = 1:nDataPoints
            xall(iDataPoint) = currMonthDepths(rowIdxs(iDataPoint));
            gall(iDataPoint) = currMonthDataType(rowIdxs(iDataPoint));
        end
        
        classicRawProfileValues(1:nDataPoints,iMonth,iLoc) = yall.*MOLAR_MASS_CARBON; % mmol C m-2 d-1 --> mg C m-2 d-1
        classicRawProfileErrSys(1:nDataPoints,iMonth,iLoc) = currMonthErrSys(1:nDataPoints).*MOLAR_MASS_CARBON; % mmol C m-2 d-1 --> mg C m-2 d-1
        classicRawProfileErrRand(1:nDataPoints,iMonth,iLoc) = currMonthErrRand(1:nDataPoints).*MOLAR_MASS_CARBON; % mmol C m-2 d-1 --> mg C m-2 d-1
        classicRawProfileDepths(1:nDataPoints,iMonth,iLoc) = xall;
        classicRawProfileDataType(1:nDataPoints,iMonth,iLoc) = gall;

    end % iMonth 
end % iLoc

% Average data per unique depth and propagate error.
 
classicMonthlyProfileAvg    = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,nLocs); % mean
classicMonthlyProfileN      = zeros(MAX_NUM_DEPTHS_PER_PROFILE,12,nLocs); % number of values
classicMonthlyProfileErrTot = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,nLocs); % total error
classicMonthlyProfileDepths = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,nLocs);  
nUniqueObsDepths            = NaN(12,nLocs); 

for iLoc = 1:nLocs
    for iMonth = 1:12

        % Get values
        valso = squeeze(classicRawProfileValues(:,iMonth,iLoc));
        errssyso = squeeze(classicRawProfileErrSys(:,iMonth,iLoc));
        errsrando = squeeze(classicRawProfileErrRand(:,iMonth,iLoc));
        zo = squeeze(classicRawProfileDepths(:,iMonth,iLoc));
        zo(zo==0) = NaN;
        uniqueObsDepths = unique(zo);
        nUniqueObsDepths(iMonth,iLoc) = sum(~isnan(uniqueObsDepths));
        
        % Calculate average value and error per unique depth
        for iUniqueDepth = 1:nUniqueObsDepths(iMonth,iLoc)
            idxThisUniqueDepth = zo == uniqueObsDepths(iUniqueDepth);
            fluxesInThisUniqueDepth = valso(idxThisUniqueDepth);
            sysErrsInThisUniqueDepth = errssyso(idxThisUniqueDepth);
            randErrsInThisUniqueDepth = errsrando(idxThisUniqueDepth);
            
            % Depth value
            classicMonthlyProfileDepths(iUniqueDepth,iMonth,iLoc) = uniqueObsDepths(iUniqueDepth);
            
            % Calculate the average and get the number of associated data points
            classicMonthlyProfileAvg(iUniqueDepth,iMonth,iLoc) = mean(fluxesInThisUniqueDepth);
            classicMonthlyProfileN(iUniqueDepth,iMonth,iLoc) = sum(fluxesInThisUniqueDepth >= 0,'omitnan');
            
            % Type A uncertainty: two options,
            % (a) if no errors associated to the samples, calculate as standard deviation/sqrt(N)
            if (sum(randErrsInThisUniqueDepth) == 0)
                typeAuncertainty = std(fluxesInThisUniqueDepth)./sqrt(classicMonthlyProfileN(iUniqueDepth,iMonth,iLoc)); 
            % (b) if error associated to the samples, propagate error
            else
                [~,~,~,f_MID,f_UB,~,~] = worstcase(@(x) mean(x),fluxesInThisUniqueDepth,randErrsInThisUniqueDepth);
                typeAuncertainty = f_UB - f_MID;
            end

            % Type B uncertainty, or systematic error 
            % Error propagation of the systematic errors
            [~,~,~,f_MID,f_UB,~,~] = worstcase(@(x) mean(x),fluxesInThisUniqueDepth,sysErrsInThisUniqueDepth);
            typeBuncertainty = f_UB - f_MID;

            % Net error: sum of type A and type B errors in the quadrature
            classicMonthlyProfileErrTot(iUniqueDepth,iMonth,iLoc) =... 
                sqrt(typeAuncertainty^2 + typeBuncertainty^2);
  
        end % iUniqueDepth
    end % iMonth
end % iLoc

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 6 - BIN DATA ANNUALLY BY DEPTH HORIZON AND PROPAGATE ERROR
% -------------------------------------------------------------------------

classicAnnualDhAvg    = NaN(3,nLocs); % weighted mean
classicAnnualDhErrTot = NaN(3,nLocs);
classicAnnualDhN      = NaN(3,nLocs);
classicAnnualDhMin    = NaN(3,nLocs);
classicAnnualDhMax    = NaN(3,nLocs);

for iLoc = 1:nLocs
    for iDh = 1:3
        nSamples = squeeze(classicMonthlyDhN(iDh,:,iLoc)); % nz x 12 x nLocs
        
        if any(nSamples)
            classicAnnualDhN(iDh,iLoc) = sum(nSamples);
            
            vals = squeeze(classicMonthlyDhAvg(iDh,:,iLoc));
            errTot = squeeze(classicMonthlyDhErrTot(iDh,:,iLoc));

            % Set to 0 positions with no samples
            vals(nSamples == 0) = 0;
            errTot(nSamples == 0) = 0;

            % Calculate weighted mean (mw = ((mA*nA)+(mB*nB)+(mC*nC))/(nA+nB+nC))
            paramsWeightedAverage = vals.*nSamples;  
            calculateWeightedAvg = @(x) sum(x) ./ sum(nSamples(:), 'omitnan');
            classicAnnualDhAvg(iDh,iLoc) = calculateWeightedAvg(paramsWeightedAverage);

            % Error propagation using worstcase with the function handle and the parameters
            [~,~,~,f_MID,f_UB,~,~] = ...
                worstcase(@(paramsWeightedAverage) calculateWeightedAvg(paramsWeightedAverage),...
                paramsWeightedAverage', errTot');
            classicAnnualDhErrTot(iDh,iLoc) = f_UB - f_MID;
            
            % Min and max values
            classicAnnualDhMin(iDh,iLoc) = min(classicMonthlyDhAvg(iDh,:,iLoc));
            classicAnnualDhMax(iDh,iLoc) = max(classicMonthlyDhAvg(iDh,:,iLoc));
        end
    end
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 7 - BIN DATA ANNUALLY BY UNIQUE DEPTH AND PROPAGATE ERROR
% -------------------------------------------------------------------------

classicAnnualProfileAvg    = NaN(MAX_NUM_DEPTHS_PER_PROFILE,nLocs); % weighted mean
classicAnnualProfileErrTot = NaN(MAX_NUM_DEPTHS_PER_PROFILE,nLocs);
classicAnnualProfileN      = NaN(MAX_NUM_DEPTHS_PER_PROFILE,nLocs);
classicAnnualProfileMin    = NaN(MAX_NUM_DEPTHS_PER_PROFILE,nLocs);
classicAnnualProfileMax    = NaN(MAX_NUM_DEPTHS_PER_PROFILE,nLocs);
classicAnnualProfileDepths = NaN(MAX_NUM_DEPTHS_PER_PROFILE,nLocs);  

for iLoc = 1:nLocs

    % Get values
    vals = squeeze(classicMonthlyProfileAvg(:,:,iLoc));
    errTot = squeeze(classicMonthlyProfileErrTot(:,:,iLoc));
    nSamples = squeeze(classicMonthlyProfileN(:,:,iLoc));

    zo = squeeze(classicMonthlyProfileDepths(:,:,iLoc));
    zo(zo==0) = NaN;
    uniqueObsDepths = unique(zo);
    nUniqueObsDepths = sum(~isnan(uniqueObsDepths));

    % Calculate average value and error per unique depth
    for iUniqueDepth = 1:nUniqueObsDepths

        idxThisUniqueDepth = zo == uniqueObsDepths(iUniqueDepth);
        fluxesInThisUniqueDepth = vals(idxThisUniqueDepth);
        errsInThisUniqueDepth = errTot(idxThisUniqueDepth);
        nSamplesInThisUniqueDepth = nSamples(idxThisUniqueDepth);
        classicAnnualProfileN(iUniqueDepth,iLoc) = sum(nSamplesInThisUniqueDepth);

        % Depth value
        classicAnnualProfileDepths(iUniqueDepth,iLoc) = uniqueObsDepths(iUniqueDepth);

        % Calculate weighted mean (mw = ((mA*nA)+(mB*nB)+(mC*nC))/(nA+nB+nC))
        paramsWeightedAverage = fluxesInThisUniqueDepth.*nSamplesInThisUniqueDepth;  
        calculateWeightedAvg = @(x) sum(x) ./ sum(nSamplesInThisUniqueDepth(:), 'omitnan');
        classicAnnualProfileAvg(iUniqueDepth,iLoc) = calculateWeightedAvg(paramsWeightedAverage);

        % Error propagation using worstcase with the function handle and the parameters
        [~,~,~,f_MID,f_UB,~,~] = ...
            worstcase(@(paramsWeightedAverage) calculateWeightedAvg(paramsWeightedAverage),...
            paramsWeightedAverage, errsInThisUniqueDepth);
        classicAnnualProfileErrTot(iUniqueDepth,iLoc) = f_UB - f_MID;

        % Min and max values
        classicAnnualProfileMin(iUniqueDepth,iLoc) = min(fluxesInThisUniqueDepth);
        classicAnnualProfileMax(iUniqueDepth,iLoc) = max(fluxesInThisUniqueDepth);

    end % iUniqueDepth
end % iLoc

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 8 - CALCULATE THE NUMBER OF DATA POINTS BASED ON VARIOUS CRITERIA
% -------------------------------------------------------------------------

% Entries for depth horizon analysis
classicMonthlyDhNbyMethod = zeros(12,nLocs,3,2); % same as classicMonthlyDhN but split by collection method
for iLoc = 1:nLocs
    for iMonth = 1:12
        for iDh = 1:3
            myvals = classicRawDhValues_cell{iDh,iMonth,iLoc};
            mytypes = classicRawDhDataType_cell{iDh,iMonth,iLoc};
            if ~isempty(myvals)
                thetypes = textscan(mytypes,'%s');
                unnestthetypes = [thetypes{:}];
                ntrap = nnz(strcmp(unnestthetypes,'trap'));
                nradio = nnz(strcmp(unnestthetypes,'radionuclide'));
                classicMonthlyDhNbyMethod(iMonth,iLoc,iDh,1) = ntrap;
                classicMonthlyDhNbyMethod(iMonth,iLoc,iDh,2) = nradio;
            end
        end
    end
end

nObsForDhAnalysis = sum(squeeze(classicMonthlyDhNbyMethod(:,:,:,:)),'all','omitnan'); % same as sum(squeeze(dataMonthlyN(:,:,:)),'all','omitnan');
fprintf('\nWe have %d data points FOR SUMMARIES BY DEPTH HORIZON, of which', nObsForDhAnalysis)

fracBatsDataPointsDhAnalysis  = (sum(squeeze(classicMonthlyDhNbyMethod(:,4,:,:)),'all','omitnan')/nObsForDhAnalysis)*100;
fracOspDataPointsDhAnalysis   = (sum(squeeze(classicMonthlyDhNbyMethod(:,2,:,:)),'all','omitnan')/nObsForDhAnalysis)*100;
fracPapsoDataPointsDhAnalysis = (sum(squeeze(classicMonthlyDhNbyMethod(:,3,:,:)),'all','omitnan')/nObsForDhAnalysis)*100;
fracEqpacDataPointsDhAnalysis = (sum(squeeze(classicMonthlyDhNbyMethod(:,1,:,:)),'all','omitnan')/nObsForDhAnalysis)*100;
fracHotalohaDataPointsDhAnalysis = (sum(squeeze(classicMonthlyDhNbyMethod(:,5,:,:)),'all','omitnan')/nObsForDhAnalysis)*100;
fracHausgartenDataPointsDhAnalysis = (sum(squeeze(classicMonthlyDhNbyMethod(:,6,:,:)),'all','omitnan')/nObsForDhAnalysis)*100;

fprintf('\n%0.1f%% data points at BATS/OFP,', fracBatsDataPointsDhAnalysis)
fprintf('\n%0.1f%% data points at HOT/ALOHA, and', fracHotalohaDataPointsDhAnalysis)
fprintf('\n%0.1f%% data points at OSP,', fracOspDataPointsDhAnalysis)
fprintf('\n%0.1f%% data points at PAP-SO,', fracPapsoDataPointsDhAnalysis)
fprintf('\n%0.1f%% data points at HAUSGARTEN.', fracHausgartenDataPointsDhAnalysis)
fprintf('\n%0.1f%% data points at EqPac,', fracEqpacDataPointsDhAnalysis)

nRadionuclideDataPoints = sum(squeeze(classicMonthlyDhNbyMethod(:,:,:,2)),'all','omitnan');
fracRadionuclideDataPoints = (nRadionuclideDataPoints/nObsForDhAnalysis)*100;
fprintf('\nWe have %0.1f%% radionuclide data points for summaries by depth horizon.', fracRadionuclideDataPoints)

% Compare with length of the raw data set
nEntriesRawDatasetByLoc = zeros(nLocs,1);
for iLoc = 1:nLocs
    currStatData = D.POC_mmol_m2_d(D.tag == STATION_NAMES{iLoc},:);
    nEntriesRawDatasetByLoc(iLoc) = nnz(~isnan(currStatData));
end
nEntriesRawDatasetTotal = sum(nEntriesRawDatasetByLoc); % = nnz(~isnan(dataRaw.POC_mmol_m2_d));
fprintf('\nWe have %d data points IN TOTAL, of which', nEntriesRawDatasetTotal)

fracBatsDataPointsRaw   = (nEntriesRawDatasetByLoc(4)/nEntriesRawDatasetTotal)*100;
fracOspDataPointsRaw    = (nEntriesRawDatasetByLoc(2)/nEntriesRawDatasetTotal)*100;
fracPapsoDataPointsRaw  = (nEntriesRawDatasetByLoc(3)/nEntriesRawDatasetTotal)*100;
fracEqpacDataPointsRaw  = (nEntriesRawDatasetByLoc(1)/nEntriesRawDatasetTotal)*100;
fracHotalohaDataPointsRaw  = (nEntriesRawDatasetByLoc(5)/nEntriesRawDatasetTotal)*100;
fracHusgartenDataPointsRaw  = (nEntriesRawDatasetByLoc(6)/nEntriesRawDatasetTotal)*100;

fprintf('\n%0.1f%% data points at BATS/OFP', fracBatsDataPointsRaw)
fprintf('\n%0.1f%% data points at HOT/ALOHA', fracHotalohaDataPointsRaw)
fprintf('\n%0.1f%% data points at OSP', fracOspDataPointsRaw)
fprintf('\n%0.1f%% data points at PAP-SO', fracPapsoDataPointsRaw)
fprintf('\n%0.1f%% data points at HAUSGARTEN', fracHusgartenDataPointsRaw)
fprintf('\n%0.1f%% data points at EqPac', fracEqpacDataPointsRaw)

nEntriesRawDatasetRadionuclideTotal = length(D.POC_mmol_m2_d(strcmp(D.method,'radionuclide'),:));
fracRadionuclidesDatasetRaw = (nEntriesRawDatasetRadionuclideTotal/nEntriesRawDatasetTotal)*100;
fprintf('\nWe have %0.1f%% radionuclide data points in total.\n', fracRadionuclidesDatasetRaw)

% Print number of data points for reference
tablendp = zeros(3,nLocs,2);
for iLoc = 1:nLocs
    for iDh = 1:3
        for iMethod = 1:2
            ndp = sum(squeeze(classicMonthlyDhNbyMethod(:,iLoc,iDh,iMethod)),'omitnan');
            tablendp(iDh,iLoc,iMethod) = ndp;
        end
    end
end
checkNumEntries = tablendp(1,6,2); % see in 6th loc, zeu for method 2 (radionuclide)

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 9 - SAVE THE DATA
% -------------------------------------------------------------------------

% Save monthly fluxes
save(fullfile('.','data','processed',filenameOutputFluxCompilation),...
    'classicRawProfileValues_cell','classicRawProfileErrRand_cell','classicRawProfileErrSys_cell',...
    'classicRawProfileDepths_cell','classicRawProfileDataType_cell','classicRawProfileValues',...
    'classicRawProfileErrRand','classicRawProfileErrSys','classicRawProfileDepths','classicRawProfileDataType',...
    'classicRawDhValues_cell','classicRawDhDepths_cell','classicRawDhTag_cell','classicRawDhDataType_cell',...
    'classicMonthlyProfileAvg','classicMonthlyProfileErrTot','classicMonthlyProfileN','classicMonthlyProfileDepths',...
    'classicMonthlyDhAvg','classicMonthlyDhErrTot','classicMonthlyDhN',...
    'classicAnnualDhAvg','classicAnnualDhErrTot','classicAnnualDhN',...
    'classicAnnualProfileAvg','classicAnnualProfileErrTot','classicAnnualProfileN','classicAnnualProfileDepths',...
    'nUniqueObsDepths','TRAPRAD_TABLE','-v7.3')

fprintf('\nThe POC flux compilation data have been saved correctly.\n')
