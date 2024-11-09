
% ======================================================================= %
%                                                                         %
% This script calculates the BCP mesopelagic transfer efficiency metrics  % 
% shown in Figure 4 in the paper (Martin's b coefficient, z* coefficient  % 
% and Teff). The metrics originate from five sources:                     % 
%     (i) measurements of POC flux from sediment traps and radionuclides  %
%         (compilation made for this study);                              %
%    (ii) estimations of POC flux from UVP5 data (from the compilation of % 
%         Kiko2022);                                                      %
%   (iii) published BCP metrics calculated from canonical fits to POC     %
%         flux data from sediment traps, radionuclides, UVP5 and in situ  %
%         filtration systems (Francois2002, Buesseler2009, Lam2011,       %
%         Guidi2015 and Mouw2016b);                                       % 
%    (iv) published BCP metrics calculated from statistical fits to       %
%         physical and biogeochemical data (Henson2012 and Marsay2015);   % 
%         and                                                             %
%     (v) published BCP metrics estimated using a diagnostic model        %
%         constrained by biogeochemical data (Weber2016).                 %
%                                                                         %
% This script uses these external functions:                              %
%   calculateBcpMetricsFromTrapAndRadionuclide.m - custom function        %
%   calculateBcpMetricsFromUvp.m                 - custom function        %
%   calculateBcpMetricsFromHenson2012.m          - custom function        %
%   calculateBcpMetricsFromMarsay2015.m          - custom function        %
%   generateMCparameters.m                       - from FileExchange      %
%   propagateErrorWithMC.m                       - from FileExchange      %
%   swtest.m                                     - from FileExchange      %
%                                                                         %
% The script has 8 sections:                                              %
%   Section 1 - Presets.                                                  %
%   Section 2 - BCP metrics computed by fitting canonical equations to    %
%               measurements of POC flux compiled for this study.         %
%   Section 3 - BCP metrics computed by fitting canonical equations to    %
%               UVP5-derived estimates of POC flux.                       %
%   Section 4 - Published BCP metrics computed from canonical fits to POC %
%               flux data.                                                %
%   Section 5 - Published BCP metrics calculated from statistical fits    %
%               to physical and biogeochemical data.                      %
%   Section 6 - Published BCP metrics estimated using a diagnostic model  %
%               constrained by biogeochemical data.                       %
%   Section 7 - Write out the metrics data 5D array into a .csv file      %
%               (Dataset S2).                                             %
%   Section 8 - ANOVA statistical tests.                                  %
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

% Choices to be applied to sediment trap/radionuclide data and UVP5 data
% when fitting b and z* (best choices determined after running
% 'processPocFluxFits.m')
isMeansOfMeans = 1;   % 0=calculate metrics from annual POC flux values, 1=calculate metrics as means of monthly metric means 
isLogTransformed = 0; % 0=do not log-log transform depth and POC flux data before fitting, 1=do transform 
isFluxNormalised = 0; % 0=do not normalise POC flux values to value at zref, 1=do normalise
choiceZref = 2;       % reference depth, 1=closest value to 100, 2=zeu, 3=inflexion point

% Path and filename declarations
filenameInputTimeseriesInformation = 'timeseries_station_information.mat';
filenameOutputMetricsData = 'bcpmetrics_all';
filenameOutputCsvTable = 'dataset_s2_figure4.csv';

% Load station information
load(fullfile('.','data','processed',filenameInputTimeseriesInformation))

% Set the seed for reproducibility of random number generator
rng(0); 

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

% Initialise the array to store metrics data
% 3rd dimension: 1=mean, 2=1stdev upp, 3=1stdev low 
% 4th dimension: up to 4 values provided by a publication
metricsData = NaN(NUM_LOCS,nMetrics,nPublications,3,4); 

% Define functions for later use
z0 = 100;
z1 = 1000;
funcTeff_bDep       = @(x) (z1/z0).^(-x); % x is Martin's b
funcZstar_bDep      = @(x) -(z1-z0)./log((z1/z0).^(-x)); % x is Martin's b
funcMartinb_teffDep = @(x) -(log(x)./(log(z1)-log(z0))); % x is Teff
funcZstar_teffDep   = @(x) -(z1-z0)./log(x); % x is Teff

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - BCP METRICS COMPUTED BY FITTING CANONICAL EQUATIONS TO
% OBSERVATIONS OF POC FLUX COMPILED FOR THIS STUDY 
% -------------------------------------------------------------------------

fprintf('\nInitiating calculations of BCP metrics from our compilation of POC flux...\n')
[classic] = calculateBcpMetricsFromTrapAndRadionuclide(isMeansOfMeans,...
    isLogTransformed,isFluxNormalised,choiceZref);
fprintf('\n...done.\n')
    
metricsData(:,iMartinb,iTimeSeries,1,1) = classic.martinb.ave(:);
metricsData(:,iMartinb,iTimeSeries,2,1) = classic.martinb.stdevupp(:);
metricsData(:,iMartinb,iTimeSeries,3,1) = classic.martinb.stdevlow(:);

metricsData(:,iZstar,iTimeSeries,1,1) = classic.zstar.ave(:);
metricsData(:,iZstar,iTimeSeries,2,1) = classic.zstar.stdevupp(:);
metricsData(:,iZstar,iTimeSeries,3,1) = classic.zstar.stdevlow(:);

metricsData(:,iTeff,iTimeSeries,1,1) = classic.teff.ave(:);
metricsData(:,iTeff,iTimeSeries,2,1) = classic.teff.stdevupp(:);
metricsData(:,iTeff,iTimeSeries,3,1) = classic.teff.stdevlow(:); 

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - BCP METRICS COMPUTED BY FITTING CANONICAL EQUATIONS TO 
% UVP5-DERIVED ESTIMATES OF POC FLUX
% -------------------------------------------------------------------------

fprintf('\nInitiating calculations of BCP metrics from the UVP5-derived estimates of POC flux...\n')
[uvp] = calculateBcpMetricsFromUvp(isMeansOfMeans,isLogTransformed,...
    isFluxNormalised,choiceZref);
fprintf('\n...done.\n')
    
metricsData(:,iMartinb,iUvp,1,1) = uvp.martinb.ave(:);
metricsData(:,iMartinb,iUvp,2,1) = uvp.martinb.stdevupp(:);
metricsData(:,iMartinb,iUvp,3,1) = uvp.martinb.stdevlow(:);

metricsData(:,iZstar,iUvp,1,1) = uvp.zstar.ave(:);
metricsData(:,iZstar,iUvp,2,1) = uvp.zstar.stdevupp(:);
metricsData(:,iZstar,iUvp,3,1) = uvp.zstar.stdevlow(:);

metricsData(:,iTeff,iUvp,1,1) = uvp.teff.ave(:);
metricsData(:,iTeff,iUvp,2,1) = uvp.teff.stdevupp(:);
metricsData(:,iTeff,iUvp,3,1) = uvp.teff.stdevlow(:);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - PUBLISHED BCP METRICS COMPUTED FROM CANONICAL FITS TO POC
% FLUX DATA
% -------------------------------------------------------------------------

% The following two steps are common to all publications:
% 1. Transform all types of reported uncertainties to standard deviation.
% 2. Propagate the error using Monte Carlo.

fprintf('\nInitiating calculations of BCP metrics from published data of BCP metrics...\n')

% .........................................................................

% Francois et al. (2002), Table 1, mean b

metricsData(iE,iMartinb,iF2002,1,1) = 0.66;   % EqPac-0, from data reported in Honjo et al. 1995 - EP=13.1, z=2284 m, F=1.64, z0=100
metricsData(iE,iMartinb,iF2002,1,2) = 0.81;   % ""                                               - EP=20.8, z=2284 m, F=1.64, z0=100
metricsData(iE,iMartinb,iF2002,1,3) = 0.59;   % ""                                               - EP=13.1, z=3618 m, F=1.60, z0=100
metricsData(iE,iMartinb,iF2002,1,4) = 0.71;   % ""                                               - EP=20.8, z=3618 m, F=1.60, z0=100
metricsData(iO,iMartinb,iF2002,1,1) = 0.80;   % P, from data reported in Wong et al. 1999        - EP=22.5, z=3800, F=1.23, z0=100
metricsData(iB,iMartinb,iF2002,1,1) = 0.97;   % OFP, from data reported in Conte et al. 2001     - EP=17.0, z=3200, F=0.59, z0=100
metricsData(iHa,iMartinb,iF2002,1,1) = 1.61;  % Fram Strait, from data reported in Honjo 1990

metricsData(iE,iTeff,iF2002,1,1) = funcTeff_bDep(metricsData(iE,iMartinb,iF2002,1,1));  
metricsData(iE,iTeff,iF2002,1,2) = funcTeff_bDep(metricsData(iE,iMartinb,iF2002,1,2));  
metricsData(iE,iTeff,iF2002,1,3) = funcTeff_bDep(metricsData(iE,iMartinb,iF2002,1,3));  
metricsData(iE,iTeff,iF2002,1,4) = funcTeff_bDep(metricsData(iE,iMartinb,iF2002,1,4));  
metricsData(iO,iTeff,iF2002,1,1) = funcTeff_bDep(metricsData(iO,iMartinb,iF2002,1,1)); 
metricsData(iB,iTeff,iF2002,1,1) = funcTeff_bDep(metricsData(iB,iMartinb,iF2002,1,1));
metricsData(iHa,iTeff,iF2002,1,1) = funcTeff_bDep(metricsData(iHa,iMartinb,iF2002,1,1));

metricsData(iE,iZstar,iF2002,1,1) = funcZstar_bDep(metricsData(iE,iMartinb,iF2002,1,1));
metricsData(iE,iZstar,iF2002,1,2) = funcZstar_bDep(metricsData(iE,iMartinb,iF2002,1,2));
metricsData(iE,iZstar,iF2002,1,3) = funcZstar_bDep(metricsData(iE,iMartinb,iF2002,1,3));
metricsData(iE,iZstar,iF2002,1,4) = funcZstar_bDep(metricsData(iE,iMartinb,iF2002,1,4));
metricsData(iO,iZstar,iF2002,1,1) = funcZstar_bDep(metricsData(iO,iMartinb,iF2002,1,1));
metricsData(iB,iZstar,iF2002,1,1) = funcZstar_bDep(metricsData(iB,iMartinb,iF2002,1,1));
metricsData(iHa,iZstar,iF2002,1,1) = funcZstar_bDep(metricsData(iHa,iMartinb,iF2002,1,1));

% .........................................................................

% Buesseler & Boyd (2009), Table 6, mean and standard error (SE) of b and z*

b2009.martinb.ave = [1.03, 1.05, 1.16, 1.17]; % Martin's b mean
b2009.martinb.se  = [0.08, 0.07, 0.14, 0.20]; % standard error = standard deviation/sqrt(N)
b2009.n           = [   2,    1,    1,    4]; % number of samples (reported in their Table 1)
b2009.martinb.stdev = b2009.martinb.se.*sqrt(b2009.n); % standard deviation

metricsData(iE,iMartinb,iB2009,1,1) = b2009.martinb.ave(1);                          % EQPAC, spring and fall avg. 1992, uses data from Bacon et al. (1996)
metricsData(iE,iMartinb,iB2009,2,1) = b2009.martinb.ave(1) + b2009.martinb.stdev(1); % EQPAC
metricsData(iE,iMartinb,iB2009,3,1) = b2009.martinb.ave(1) - b2009.martinb.stdev(1); % EQPAC
metricsData(iO,iMartinb,iB2009,1,1) = b2009.martinb.ave(2);                          % OSP-May, uses data from Charette et al. (1999)
metricsData(iO,iMartinb,iB2009,2,1) = b2009.martinb.ave(2) + b2009.martinb.stdev(2); % OSP-May
metricsData(iO,iMartinb,iB2009,3,1) = b2009.martinb.ave(2) - b2009.martinb.stdev(2); % OSP-May
metricsData(iO,iMartinb,iB2009,1,2) = b2009.martinb.ave(3);                          % OSP-Aug, ""
metricsData(iO,iMartinb,iB2009,2,2) = b2009.martinb.ave(3) + b2009.martinb.stdev(3); % OSP-Aug
metricsData(iO,iMartinb,iB2009,3,2) = b2009.martinb.ave(3) - b2009.martinb.stdev(3); % OSP-Aug
metricsData(iHo,iMartinb,iB2009,1,1) = b2009.martinb.ave(4);                          % HOT, Jun/Jul 2004, uses data from Buesseler et al. (2007, 2008)
metricsData(iHo,iMartinb,iB2009,2,1) = b2009.martinb.ave(4) + b2009.martinb.stdev(4); % HOT, Jun/Jul 2004, uses data from Buesseler et al. (2007, 2008)
metricsData(iHo,iMartinb,iB2009,3,1) = b2009.martinb.ave(4) - b2009.martinb.stdev(4); % HOT, Jun/Jul 2004, uses data from Buesseler et al. (2007, 2008)

b2009.zstar.ave   = [217, 78, 77, 216]; 
b2009.zstar.se    = [ 24, 20, 20,  74]; 
b2009.zstar.stdev = b2009.zstar.se.*sqrt(b2009.n); 

metricsData(iE,iZstar,iB2009,1,1) = b2009.zstar.ave(1);            
metricsData(iE,iZstar,iB2009,2,1) = b2009.zstar.ave(1) + b2009.zstar.stdev(1);    
metricsData(iE,iZstar,iB2009,3,1) = b2009.zstar.ave(1) - b2009.zstar.stdev(1);    
metricsData(iO,iZstar,iB2009,1,1) = b2009.zstar.ave(2);            
metricsData(iO,iZstar,iB2009,2,1) = b2009.zstar.ave(2) + b2009.zstar.stdev(2);   
metricsData(iO,iZstar,iB2009,3,1) = b2009.zstar.ave(2) - b2009.zstar.stdev(2);    
metricsData(iO,iZstar,iB2009,1,2) = b2009.zstar.ave(3);            
metricsData(iO,iZstar,iB2009,2,2) = b2009.zstar.ave(3) + b2009.zstar.stdev(3);    
metricsData(iO,iZstar,iB2009,3,2) = b2009.zstar.ave(3) - b2009.zstar.stdev(3); 
metricsData(iHo,iZstar,iB2009,1,1) = b2009.zstar.ave(4);            
metricsData(iHo,iZstar,iB2009,2,1) = b2009.zstar.ave(4) + b2009.zstar.stdev(4);    
metricsData(iHo,iZstar,iB2009,3,1) = b2009.zstar.ave(4) - b2009.zstar.stdev(4); 

ci = NaN(2,length(b2009.martinb.ave));
midval = NaN(1,length(b2009.martinb.ave));
for i = 1:length(b2009.martinb.ave)
    % Error propagation using MC, assume gaussian distrib of errors
    A = generateMCparameters('gaussian',[b2009.martinb.ave(i),b2009.martinb.stdev(i)],'plot',false);
    [midval(i),ci(:,i),~] = propagateErrorWithMC(funcTeff_bDep,A,'plot',false);   
end

metricsData(iE,iTeff,iB2009,1,1) = midval(1); 
metricsData(iE,iTeff,iB2009,2,1) = ci(2,1); 
metricsData(iE,iTeff,iB2009,3,1) = ci(1,1); 
metricsData(iO,iTeff,iB2009,1,1) = midval(2); 
metricsData(iO,iTeff,iB2009,2,1) = ci(2,2); 
metricsData(iO,iTeff,iB2009,3,1) = ci(1,2);
metricsData(iO,iTeff,iB2009,1,2) = midval(3); 
metricsData(iO,iTeff,iB2009,2,2) = ci(2,3); 
metricsData(iO,iTeff,iB2009,3,2) = ci(1,3); 
metricsData(iHo,iTeff,iB2009,1,1) = midval(4); 
metricsData(iHo,iTeff,iB2009,2,1) = ci(2,4); 
metricsData(iHo,iTeff,iB2009,3,1) = ci(1,4); 

% .........................................................................

% Lam et al. (2011), Suppl., mean and 95% confidence intervals

% We use the following protocol highlighted in the Cochrane Handbook to
% calculate standard deviation from confidence intervals. We assume N = 1.
% https://handbook-5-1.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm
% CI(95%): zscore = 1.96
% CI(upp) = mean + zscore * std/sqrt(N) --> mean = CI(upp) - zscore * std/sqrt(N)
% CI(low) = mean - zscore * std/sqrt(N) --> mean = CI(low) + zscore * std/sqrt(N)
% thus, 
%   CI(upp) - zscore * std/sqrt(N) = CI(low) + zscore * std/sqrt(N)
%   CI(upp) - CI(low) = 2 * (zscore * std/sqrt(N)) = 3.92 * std/sqrt(N)
%   std = (CI(upp) - CI(low)) * sqrt(N) / 3.92

l2011.martinb.ci95.upp = [1.35,  1.69, 1.98, 1.17, 1.29, 2.07, 2.45, 1.64]; % Martin's b upper 95% CI
l2011.martinb.ci95.low = [-0.70, 1.00, 1.09, 0.72, 0.41, 0.91, 0.59, 0.48]; % lower 95% CI – corrected 2nd value from 0.10 to 1.00
l2011.martinb.ave      = (l2011.martinb.ci95.upp + l2011.martinb.ci95.low)./2; % recalculate the mean (compare with [0.32,  1.34, 1.54, 0.94, 0.85])
l2011.martinb.stdev    = (l2011.martinb.ci95.upp + l2011.martinb.ci95.low)./3.92; % standard deviation (=68% CI)

metricsData(iE,iMartinb,iL2011,1,1) = l2011.martinb.ave(1);                          % EqPac EQ1, MULVFS, data from Bishop (1992) and Bishop (1999)
metricsData(iE,iMartinb,iL2011,2,1) = l2011.martinb.ave(1) + l2011.martinb.stdev(1); % EqPac EQ1
metricsData(iE,iMartinb,iL2011,3,1) = l2011.martinb.ave(1) - l2011.martinb.stdev(1); % EqPac EQ1
metricsData(iE,iMartinb,iL2011,1,2) = l2011.martinb.ave(2);                          % EqPac EQ2, MULVFS, ""
metricsData(iE,iMartinb,iL2011,2,2) = l2011.martinb.ave(2) + l2011.martinb.stdev(2); % EqPac EQ2
metricsData(iE,iMartinb,iL2011,3,2) = l2011.martinb.ave(2) - l2011.martinb.stdev(2); % EqPac EQ2
metricsData(iO,iMartinb,iL2011,1,1) = l2011.martinb.ave(3);                          % OSP Feb 96, MULVFS, data from Bishop et al. (1999)
metricsData(iO,iMartinb,iL2011,2,1) = l2011.martinb.ave(3) + l2011.martinb.stdev(3); % OSP Feb 96
metricsData(iO,iMartinb,iL2011,3,1) = l2011.martinb.ave(3) - l2011.martinb.stdev(3); % OSP Feb 96
metricsData(iO,iMartinb,iL2011,1,2) = l2011.martinb.ave(4);                          % OSP May 96, MULVFS, ""
metricsData(iO,iMartinb,iL2011,2,2) = l2011.martinb.ave(4) + l2011.martinb.stdev(4); % OSP May 96
metricsData(iO,iMartinb,iL2011,3,2) = l2011.martinb.ave(4) - l2011.martinb.stdev(4); % OSP May 96
metricsData(iO,iMartinb,iL2011,1,3) = l2011.martinb.ave(5);                          % OSP Aug 96, MULVFS, ""
metricsData(iO,iMartinb,iL2011,2,3) = l2011.martinb.ave(5) + l2011.martinb.stdev(5); % OSP Aug 96
metricsData(iO,iMartinb,iL2011,3,3) = l2011.martinb.ave(5) - l2011.martinb.stdev(5); % OSP Aug 96
metricsData(iHo,iMartinb,iL2011,1,1) = l2011.martinb.ave(6);                          % ALOHA Jun/Jul 2004, MULVFS, data from Bishop & Wood (2008)
metricsData(iHo,iMartinb,iL2011,2,1) = l2011.martinb.ave(6) + l2011.martinb.stdev(6); % ""
metricsData(iHo,iMartinb,iL2011,3,1) = l2011.martinb.ave(6) - l2011.martinb.stdev(6); % ""
metricsData(iHo,iMartinb,iL2011,1,2) = l2011.martinb.ave(7);                          % ""
metricsData(iHo,iMartinb,iL2011,2,2) = l2011.martinb.ave(7) + l2011.martinb.stdev(7); % ""
metricsData(iHo,iMartinb,iL2011,3,2) = l2011.martinb.ave(7) - l2011.martinb.stdev(7); % ""
metricsData(iHo,iMartinb,iL2011,1,3) = l2011.martinb.ave(8);                          % ""
metricsData(iHo,iMartinb,iL2011,2,3) = l2011.martinb.ave(8) + l2011.martinb.stdev(8); % ""
metricsData(iHo,iMartinb,iL2011,3,3) = l2011.martinb.ave(8) - l2011.martinb.stdev(8); % ""

ci = NaN(2,length(l2011.martinb.ave));
midval = NaN(1,length(l2011.martinb.ave));
for i = 1:length(l2011.martinb.ave)
    % Error propagation using MC, assume gaussian distrib of errors
    A = generateMCparameters('gaussian',[l2011.martinb.ave(i),l2011.martinb.stdev(i)],'plot',false);
    [midval(i),ci(:,i),funvals] = propagateErrorWithMC(funcTeff_bDep,A,'plot',false); 
end

metricsData(iE,iTeff,iL2011,1,1) = midval(1);  
metricsData(iE,iTeff,iL2011,2,1) = ci(2,1);  
metricsData(iE,iTeff,iL2011,3,1) = ci(1,1); 
metricsData(iE,iTeff,iL2011,1,2) = midval(2);
metricsData(iE,iTeff,iL2011,2,2) = ci(2,2);
metricsData(iE,iTeff,iL2011,3,2) = ci(1,2);
metricsData(iO,iTeff,iL2011,1,1) = midval(3);
metricsData(iO,iTeff,iL2011,2,1) = ci(2,3);
metricsData(iO,iTeff,iL2011,3,1) = ci(1,3);
metricsData(iO,iTeff,iL2011,1,2) = midval(4);
metricsData(iO,iTeff,iL2011,2,2) = ci(2,4);
metricsData(iO,iTeff,iL2011,3,2) = ci(1,4);
metricsData(iO,iTeff,iL2011,1,3) = midval(5);
metricsData(iO,iTeff,iL2011,2,3) = ci(2,5); 
metricsData(iO,iTeff,iL2011,3,3) = ci(1,5);
metricsData(iHo,iTeff,iL2011,1,1) = midval(6);
metricsData(iHo,iTeff,iL2011,2,1) = ci(2,6); 
metricsData(iHo,iTeff,iL2011,3,1) = ci(1,6);
metricsData(iHo,iTeff,iL2011,1,2) = midval(7);
metricsData(iHo,iTeff,iL2011,2,2) = ci(2,7); 
metricsData(iHo,iTeff,iL2011,3,2) = ci(1,7);
metricsData(iHo,iTeff,iL2011,1,3) = midval(8);
metricsData(iHo,iTeff,iL2011,2,3) = ci(2,8); 
metricsData(iHo,iTeff,iL2011,3,3) = ci(1,8);

ci = NaN(2,length(l2011.martinb.ave));
midval = NaN(1,length(l2011.martinb.ave));
for i = 1:length(l2011.martinb.ave)
    % Error propagation using MC, assume gaussian distrib of errors
    A = generateMCparameters('gaussian',[l2011.martinb.ave(i),l2011.martinb.stdev(i)],'plot',false);
    [midval(i),ci(:,i),~] = propagateErrorWithMC(funcZstar_bDep,A,'plot',false);   
end

metricsData(iE,iZstar,iL2011,1,1) = midval(1);
metricsData(iE,iZstar,iL2011,2,1) = ci(2,1); 
metricsData(iE,iZstar,iL2011,3,1) = ci(1,1);
metricsData(iE,iZstar,iL2011,1,2) = midval(2);
metricsData(iE,iZstar,iL2011,2,2) = ci(2,2);
metricsData(iE,iZstar,iL2011,3,2) = ci(1,2);
metricsData(iO,iZstar,iL2011,1,1) = midval(3);
metricsData(iO,iZstar,iL2011,2,1) = ci(2,3);
metricsData(iO,iZstar,iL2011,3,1) = ci(1,3); 
metricsData(iO,iZstar,iL2011,1,2) = midval(4);
metricsData(iO,iZstar,iL2011,2,2) = ci(2,4);
metricsData(iO,iZstar,iL2011,3,2) = ci(1,4); 
metricsData(iO,iZstar,iL2011,1,3) = midval(5);
metricsData(iO,iZstar,iL2011,2,3) = ci(2,5);
metricsData(iO,iZstar,iL2011,3,3) = ci(1,5);
metricsData(iHo,iZstar,iL2011,1,1) = midval(6);
metricsData(iHo,iZstar,iL2011,2,1) = ci(2,6);
metricsData(iHo,iZstar,iL2011,3,1) = ci(1,6);
metricsData(iHo,iZstar,iL2011,1,2) = midval(7);
metricsData(iHo,iZstar,iL2011,2,2) = ci(2,7);
metricsData(iHo,iZstar,iL2011,3,2) = ci(1,7);
metricsData(iHo,iZstar,iL2011,1,3) = midval(8);
metricsData(iHo,iZstar,iL2011,2,3) = ci(2,8);
metricsData(iHo,iZstar,iL2011,3,3) = ci(1,8);

% .........................................................................

% Guidi et al. (2015), Suppl., median and interquartile range

% In order to estimate the sample mean and standard deviation from the 
% median and interquartile range, we've used the method proposed in 
% Wan et al. (2014), case 3, which updates the method proposed in the 
% Cochrane Handbook. This method performs well for normal and skewed data.
% Paper here: https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-14-135

g2015.martinb.med = [0.69, 0.82, 0.84, 0.95, 0.68, 0.89, 1.35]; % Martin's b median
g2015.n           = [  60,   13,   27,   25,    6,   13,    6]; % number of samples
g2015.martinb.q3  = [0.85, 1.02, 1.05, 1.07, 0.92, 1.11, 1.80]; % 75% of b
g2015.martinb.q1  = [0.59, 0.70, 0.63, 0.76, 0.39, 0.56, 0.88]; % 25% of b
g2015.q = (g2015.n-1)./4; % after the definition N = 4Q + 1
g2015.mu_n        = [1.318, 1.206, 1.284, 1.274, 0.990, 1.206, 0.990]; % tabulated (Table 2 in Wan et al. 2014, only if q <= 50)
g2015.martinb.ave = (g2015.martinb.q1 + g2015.martinb.med + g2015.martinb.q3)./3; % mean
g2015.martinb.stdev = (g2015.martinb.q3 - g2015.martinb.q1)./g2015.mu_n; % standard deviation

metricsData(iE,iMartinb,iG2015,1,1) = g2015.martinb.ave(1);                          % PEQD, from 234Th-Fth and sediment trap-Fs
metricsData(iE,iMartinb,iG2015,2,1) = g2015.martinb.ave(1) + g2015.martinb.stdev(1); % PEQD
metricsData(iE,iMartinb,iG2015,3,1) = g2015.martinb.ave(1) - g2015.martinb.stdev(1); % PEQD
metricsData(iO,iMartinb,iG2015,1,1) = g2015.martinb.ave(2);                          % PSAE, from 234Th-Fth and sediment trap-Fs
metricsData(iO,iMartinb,iG2015,2,1) = g2015.martinb.ave(2) + g2015.martinb.stdev(2); % PSAE
metricsData(iO,iMartinb,iG2015,3,1) = g2015.martinb.ave(2) - g2015.martinb.stdev(2); % PSAE
metricsData(iP,iMartinb,iG2015,1,1) = g2015.martinb.ave(3);                          % NADR, from UVP profiles
metricsData(iP,iMartinb,iG2015,2,1) = g2015.martinb.ave(3) + g2015.martinb.stdev(3); % NADR
metricsData(iP,iMartinb,iG2015,3,1) = g2015.martinb.ave(3) - g2015.martinb.stdev(3); % NADR
metricsData(iP,iMartinb,iG2015,1,2) = g2015.martinb.ave(4);                          % NADR, from 234Th-Fth and sediment trap-Fs 
metricsData(iP,iMartinb,iG2015,2,2) = g2015.martinb.ave(4) + g2015.martinb.stdev(4); % NADR
metricsData(iP,iMartinb,iG2015,3,2) = g2015.martinb.ave(4) - g2015.martinb.stdev(4); % NADR
metricsData(iB,iMartinb,iG2015,1,1) = g2015.martinb.ave(5);                          % NASW, from UVP profiles
metricsData(iB,iMartinb,iG2015,2,1) = g2015.martinb.ave(5) + g2015.martinb.stdev(5); % NASW
metricsData(iB,iMartinb,iG2015,3,1) = g2015.martinb.ave(5) - g2015.martinb.stdev(5); % NASW
metricsData(iHo,iMartinb,iG2015,1,1) = g2015.martinb.ave(6);                          % NPTG, from UVP profiles  
metricsData(iHo,iMartinb,iG2015,2,1) = g2015.martinb.ave(6) + g2015.martinb.stdev(6); % NPTG
metricsData(iHo,iMartinb,iG2015,3,1) = g2015.martinb.ave(6) - g2015.martinb.stdev(6); % NPTG
metricsData(iHo,iMartinb,iG2015,1,2) = g2015.martinb.ave(7);                          % NPTG, from 234Th-Fth and sediment trap-Fs
metricsData(iHo,iMartinb,iG2015,2,2) = g2015.martinb.ave(7) + g2015.martinb.stdev(7); % NPTG
metricsData(iHo,iMartinb,iG2015,3,2) = g2015.martinb.ave(7) - g2015.martinb.stdev(7); % NPTG

ci = NaN(2,length(g2015.martinb.ave));
midval = NaN(1,length(g2015.martinb.ave));
for i = 1:length(g2015.martinb.ave)
    % Error propagation using MC, assume gaussian distrib of errors
    A = generateMCparameters('gaussian',[g2015.martinb.ave(i),g2015.martinb.stdev(i)],'plot',false);
    [midval(i),ci(:,i),funvals] = propagateErrorWithMC(funcTeff_bDep,A,'plot',false);   
end

metricsData(iE,iTeff,iG2015,1,1) = midval(1);   
metricsData(iE,iTeff,iG2015,2,1) = ci(2,1);
metricsData(iE,iTeff,iG2015,3,1) = ci(1,1);
metricsData(iO,iTeff,iG2015,1,1) = midval(2);   
metricsData(iO,iTeff,iG2015,2,1) = ci(2,2);
metricsData(iO,iTeff,iG2015,3,1) = ci(1,2);
metricsData(iP,iTeff,iG2015,1,1) = midval(3);   
metricsData(iP,iTeff,iG2015,2,1) = ci(2,3);
metricsData(iP,iTeff,iG2015,3,1) = ci(1,3);
metricsData(iP,iTeff,iG2015,1,2) = midval(4);   
metricsData(iP,iTeff,iG2015,2,2) = ci(2,4);
metricsData(iP,iTeff,iG2015,3,2) = ci(1,4);
metricsData(iB,iTeff,iG2015,1,1) = midval(5);  
metricsData(iB,iTeff,iG2015,2,1) = ci(2,5);   
metricsData(iB,iTeff,iG2015,3,1) = ci(1,5);
metricsData(iHo,iTeff,iG2015,1,1) = midval(6);  
metricsData(iHo,iTeff,iG2015,2,1) = ci(2,6);   
metricsData(iHo,iTeff,iG2015,3,1) = ci(1,6);
metricsData(iHo,iTeff,iG2015,1,2) = midval(7);  
metricsData(iHo,iTeff,iG2015,2,2) = ci(2,7);   
metricsData(iHo,iTeff,iG2015,3,2) = ci(1,7);

ci = NaN(2,length(g2015.martinb.ave));
midval = NaN(1,length(g2015.martinb.ave));
for i = 1:length(g2015.martinb.ave)
    % Error propagation using MC, assume gaussian distrib of errors
    A = generateMCparameters('gaussian',[g2015.martinb.ave(i),g2015.martinb.stdev(i)],'plot',false);
    [midval(i),ci(:,i),~] = propagateErrorWithMC(funcZstar_bDep,A,'plot',false);   
end

metricsData(iE,iZstar,iG2015,1,1) = midval(1);
metricsData(iE,iZstar,iG2015,2,1) = ci(2,1);
metricsData(iE,iZstar,iG2015,3,1) = ci(1,1);
metricsData(iO,iZstar,iG2015,1,1) = midval(2);
metricsData(iO,iZstar,iG2015,2,1) = ci(2,2);
metricsData(iO,iZstar,iG2015,3,1) = ci(1,2);
metricsData(iP,iZstar,iG2015,1,1) = midval(3);
metricsData(iP,iZstar,iG2015,2,1) = ci(2,3);
metricsData(iP,iZstar,iG2015,3,1) = ci(1,3);
metricsData(iP,iZstar,iG2015,1,2) = midval(4);
metricsData(iP,iZstar,iG2015,2,2) = ci(2,4);
metricsData(iP,iZstar,iG2015,3,2) = ci(1,4);
metricsData(iB,iZstar,iG2015,1,1) = midval(5);
metricsData(iB,iZstar,iG2015,2,1) = ci(2,5);
metricsData(iB,iZstar,iG2015,3,1) = ci(1,5);
metricsData(iHo,iZstar,iG2015,1,1) = midval(6);
metricsData(iHo,iZstar,iG2015,2,1) = ci(2,6);
metricsData(iHo,iZstar,iG2015,3,1) = ci(1,6);
metricsData(iHo,iZstar,iG2015,1,2) = midval(7);
metricsData(iHo,iZstar,iG2015,2,2) = ci(2,7);
metricsData(iHo,iZstar,iG2015,3,2) = ci(1,7);

% .........................................................................

% Mouw et al. (2016b), Fig. 3 and 9 (digitised)

m2016.zstar.upp = [343, 51, 273, 184]; % z* mean + stdev upp
m2016.zstar.low = [269, 27, 213, 132]; % z* mean - stdev upp
m2016.zstar.ave = (m2016.zstar.upp + m2016.zstar.low)./2; % mean
m2016.zstar.stdev = m2016.zstar.upp - m2016.zstar.ave; % standard deviation

m2016.oneMinusAlpha.upp = [0.078, 0.243, 0.054, 0.097]; % (1-alpha) mean + stdev upp
m2016.oneMinusAlpha.low = [0.044, 0.095, 0.040, 0.047]; % (1-alpha) mean + stdev low
m2016.oneMinusAlpha.ave = (m2016.oneMinusAlpha.upp + m2016.oneMinusAlpha.low)./2; % mean
m2016.oneMinusAlpha.stdev = m2016.oneMinusAlpha.upp - m2016.oneMinusAlpha.ave; % standard deviation

metricsData(iE,iZstar,iM2016,1,1) = m2016.zstar.ave(1);                         % PEQD province, digit. Fig. 9, based on data from dataset in Mouw et al. (2016)
metricsData(iE,iZstar,iM2016,2,1) = m2016.zstar.ave(1) + m2016.zstar.stdev(1);  % PEQD province
metricsData(iE,iZstar,iM2016,3,1) = m2016.zstar.ave(1) - m2016.zstar.stdev(1);  % PEQD province
metricsData(iO,iZstar,iM2016,1,1) = m2016.zstar.ave(2);                         % OSP, digit. Fig. 3, time-series data
metricsData(iO,iZstar,iM2016,2,1) = m2016.zstar.ave(2) + m2016.zstar.stdev(2);  % OSP
metricsData(iO,iZstar,iM2016,3,1) = m2016.zstar.ave(2) - m2016.zstar.stdev(2);  % OSP
metricsData(iP,iZstar,iM2016,1,1) = m2016.zstar.ave(3);                         % NADR province, digit. Fig. 9
metricsData(iP,iZstar,iM2016,2,1) = m2016.zstar.ave(3) + m2016.zstar.stdev(3);  % NADR province
metricsData(iP,iZstar,iM2016,3,1) = m2016.zstar.ave(3) - m2016.zstar.stdev(3);  % NADR province
metricsData(iB,iZstar,iM2016,1,1) = m2016.zstar.ave(4);                         % BATS, digit. Fig. 3, time-series data
metricsData(iB,iZstar,iM2016,2,1) = m2016.zstar.ave(4) + m2016.zstar.stdev(4);  % BATS
metricsData(iB,iZstar,iM2016,3,1) = m2016.zstar.ave(4) - m2016.zstar.stdev(4);  % BATS

% Check Eq. 1 & 2 in Mouw et al. (2016b) for the definition of Teff (they 
% have their own equation), where Teff is the term in between square
% brackets in Eq. 1)

funcTeff_m2016    = @(x) (1-x(2)).*exp(-(z1-z0)./x(1)) + x(2); % x(1) is z* and x(2) is (1-alpha)
funcMartinb_m2016 = @(x) -(log((1-x(2)).*exp(-(z1-z0)./x(1)) + x(2))./(log(z1)-log(z0))); % x(1) is z* and x(2) is (1-alpha)

midval_teff    = NaN(1,5);
ci_teff        = NaN(2,5);
midval_martinb = NaN(1,5);
ci_martinb     = NaN(2,5);
for i = 1:length(m2016.zstar.upp)
    % Error propagation using MC, assume gaussian distrib of errors
    A = generateMCparameters('gaussian',[m2016.zstar.ave(i),m2016.zstar.stdev(i)],'plot',false);
    B = generateMCparameters('gaussian',[m2016.oneMinusAlpha.ave(i),m2016.oneMinusAlpha.stdev(i)],'plot',false);
    % We don't want 0 values, replace by the non-zero minimum 
    A(A<0) = min(A(A>0));
    B(B<0) = min(B(B>0));
    [midval_teff(i),ci_teff(:,i),~] = propagateErrorWithMC(funcTeff_m2016,[A;B],'plot',false);  
    [midval_martinb(i),ci_martinb(:,i),~] = propagateErrorWithMC(funcMartinb_m2016,[A;B],'plot',false);
end

metricsData(iE,iTeff,iM2016,1,1) = midval_teff(1);
metricsData(iE,iTeff,iM2016,2,1) = ci_teff(2,1);
metricsData(iE,iTeff,iM2016,3,1) = ci_teff(1,1);
metricsData(iO,iTeff,iM2016,1,1) = midval_teff(2);
metricsData(iO,iTeff,iM2016,2,1) = ci_teff(2,2);
metricsData(iO,iTeff,iM2016,3,1) = ci_teff(1,2);
metricsData(iP,iTeff,iM2016,1,1) = midval_teff(3);
metricsData(iP,iTeff,iM2016,2,1) = ci_teff(2,3);
metricsData(iP,iTeff,iM2016,3,1) = ci_teff(1,3);
metricsData(iB,iTeff,iM2016,1,1) = midval_teff(4);
metricsData(iB,iTeff,iM2016,2,1) = ci_teff(2,4); 
metricsData(iB,iTeff,iM2016,3,1) = ci_teff(1,4); 

metricsData(iE,iMartinb,iM2016,1,1) = midval_martinb(1);
metricsData(iE,iMartinb,iM2016,2,1) = ci_martinb(2,1);
metricsData(iE,iMartinb,iM2016,3,1) = ci_martinb(1,1);
metricsData(iO,iMartinb,iM2016,1,1) = midval_martinb(2);
metricsData(iO,iMartinb,iM2016,2,1) = ci_martinb(2,2);
metricsData(iO,iMartinb,iM2016,3,1) = ci_martinb(1,2);
metricsData(iP,iMartinb,iM2016,1,1) = midval_martinb(3);
metricsData(iP,iMartinb,iM2016,2,1) = ci_martinb(2,3);
metricsData(iP,iMartinb,iM2016,3,1) = ci_martinb(1,3);
metricsData(iB,iMartinb,iM2016,1,1) = midval_martinb(4);
metricsData(iB,iMartinb,iM2016,2,1) = ci_martinb(2,4);
metricsData(iB,iMartinb,iM2016,3,1) = ci_martinb(1,4);

fprintf('...done.\n')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - PUBLISHED BCP METRICS COMPUTED FROM STATISTICAL FITS TO
% PHYSICAL AND BIOGEOCHEMICAL DATA
% -------------------------------------------------------------------------

% Henson et al. (2012) 

fprintf('\nInitiating calculations of BCP metrics following Henson et al. (2012)...\n')
[henson2012] = calculateBcpMetricsFromHenson2012(0,LOC_LATS,LOC_LONS);
fprintf('\n...done.\n')

metricsData(:,iMartinb,iH2012,1,1) = henson2012.martinb.ave(:);
metricsData(:,iMartinb,iH2012,2,1) = henson2012.martinb.stdevupp(:);
metricsData(:,iMartinb,iH2012,3,1) = henson2012.martinb.stdevlow(:);
metricsData(:,iMartinb,iH2012,4,1) = henson2012.martinb.max(:);
metricsData(:,iMartinb,iH2012,5,1) = henson2012.martinb.min(:);

metricsData(:,iZstar,iH2012,1,1) = henson2012.zstar.ave(:);
metricsData(:,iZstar,iH2012,2,1) = henson2012.zstar.stdevupp(:);
metricsData(:,iZstar,iH2012,3,1) = henson2012.zstar.stdevlow(:);
metricsData(:,iZstar,iH2012,4,1) = henson2012.zstar.max(:);
metricsData(:,iZstar,iH2012,5,1) = henson2012.zstar.min(:);

metricsData(:,iTeff,iH2012,1,1) = henson2012.teff1000.ave(:);
metricsData(:,iTeff,iH2012,2,1) = henson2012.teff1000.stdevupp(:);
metricsData(:,iTeff,iH2012,3,1) = henson2012.teff1000.stdevlow(:);
metricsData(:,iTeff,iH2012,4,1) = henson2012.teff1000.max(:); 
metricsData(:,iTeff,iH2012,5,1) = henson2012.teff1000.min(:); 

% .........................................................................

% Marsay et al. (2015)

fprintf('\nInitiating calculations of BCP metrics following Marsay et al. (2015)...\n')
[marsay2015] = calculateBcpMetricsFromMarsay2015(0,LOC_LATS,LOC_LONS);
fprintf('\n...done.\n')

metricsData(:,iMartinb,iM2015,1,1) = marsay2015.martinb.ave(:);
metricsData(:,iMartinb,iM2015,2,1) = marsay2015.martinb.stdevupp(:);
metricsData(:,iMartinb,iM2015,3,1) = marsay2015.martinb.stdevlow(:);
metricsData(:,iMartinb,iM2015,4,1) = marsay2015.martinb.max(:);
metricsData(:,iMartinb,iM2015,5,1) = marsay2015.martinb.min(:);

metricsData(:,iZstar,iM2015,1,1) = marsay2015.zstar.ave(:);
metricsData(:,iZstar,iM2015,2,1) = marsay2015.zstar.stdevupp(:);
metricsData(:,iZstar,iM2015,3,1) = marsay2015.zstar.stdevlow(:); 
metricsData(:,iZstar,iM2015,4,1) = marsay2015.zstar.max(:); 
metricsData(:,iZstar,iM2015,5,1) = marsay2015.zstar.min(:); 

metricsData(:,iTeff,iM2015,1,1) = marsay2015.teff1000.ave(:);
metricsData(:,iTeff,iM2015,2,1) = marsay2015.teff1000.stdevupp(:);
metricsData(:,iTeff,iM2015,3,1) = marsay2015.teff1000.stdevlow(:);
metricsData(:,iTeff,iM2015,4,1) = marsay2015.teff1000.max(:); 
metricsData(:,iTeff,iM2015,5,1) = marsay2015.teff1000.min(:); 

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 6 - PUBLISHED BCP METRICS ESTIMATED USING A DIAGNOSTIC MODEL
% CONSTRAINED BY BIOGEOCHEMICAL DATA
% -------------------------------------------------------------------------

fprintf('\nInitiating calculations of BCP metrics derived from the metrics of Weber et al. (2016)...\n')

% Weber et al. (2016) outputs, "Teff zeu -> 1000 m", which I have obtained
% by digitising their Fig. 3C

w2016.teffZeuto1000.upp = [0.193, 0.167, 0.298, 0.261, 0.368, 0.349, 0.177, 0.161, 0.056, 0.055]; % Teff mean + stdv upp
% w2016.teffZeuto1000.low = [0.132, 0.127, 0.159, 0.184, 0.117, 0.161, -0.0441, -0.0407]; % Teff mean - stdv low
w2016.teffZeuto1000.low = [0.132, 0.127, 0.159, 0.184, 0.117, 0.161,     0,     0,     0,     0];
w2016.teffZeuto1000.ave = (w2016.teffZeuto1000.upp + w2016.teffZeuto1000.low)./2; % mean
w2016.teffZeuto1000.stdev = w2016.teffZeuto1000.upp - w2016.teffZeuto1000.ave; % standard deviation

metricsData(iE,iTeff,iW2016,1,1) = w2016.teffZeuto1000.ave(1);                                % ETP, satellite export
metricsData(iE,iTeff,iW2016,2,1) = w2016.teffZeuto1000.ave(1) + w2016.teffZeuto1000.stdev(1); % ETP
metricsData(iE,iTeff,iW2016,3,1) = w2016.teffZeuto1000.ave(1) - w2016.teffZeuto1000.stdev(1); % ETP
metricsData(iE,iTeff,iW2016,1,2) = w2016.teffZeuto1000.ave(2);                                % ETP, CMIP export
metricsData(iE,iTeff,iW2016,2,2) = w2016.teffZeuto1000.ave(2) + w2016.teffZeuto1000.stdev(2); % ETP
metricsData(iE,iTeff,iW2016,3,2) = w2016.teffZeuto1000.ave(2) - w2016.teffZeuto1000.stdev(2); % ETP
metricsData(iO,iTeff,iW2016,1,1) = w2016.teffZeuto1000.ave(3);                                % NP, satellite export
metricsData(iO,iTeff,iW2016,2,1) = w2016.teffZeuto1000.ave(3) + w2016.teffZeuto1000.stdev(3); % NP
metricsData(iO,iTeff,iW2016,3,1) = w2016.teffZeuto1000.ave(3) - w2016.teffZeuto1000.stdev(3); % NP
metricsData(iO,iTeff,iW2016,1,2) = w2016.teffZeuto1000.ave(4);                                % NP, CMIP export
metricsData(iO,iTeff,iW2016,2,2) = w2016.teffZeuto1000.ave(4) + w2016.teffZeuto1000.stdev(4); % NP
metricsData(iO,iTeff,iW2016,3,2) = w2016.teffZeuto1000.ave(4) - w2016.teffZeuto1000.stdev(4); % NP
metricsData(iP,iTeff,iW2016,1,1) = w2016.teffZeuto1000.ave(5);                                % NA, satellite export
metricsData(iP,iTeff,iW2016,2,1) = w2016.teffZeuto1000.ave(5) + w2016.teffZeuto1000.stdev(5); % NA
metricsData(iP,iTeff,iW2016,3,1) = w2016.teffZeuto1000.ave(5) - w2016.teffZeuto1000.stdev(5); % NA
metricsData(iP,iTeff,iW2016,1,2) = w2016.teffZeuto1000.ave(6);                                % NA, CMIP export
metricsData(iP,iTeff,iW2016,2,2) = w2016.teffZeuto1000.ave(6) + w2016.teffZeuto1000.stdev(6); % NA
metricsData(iP,iTeff,iW2016,3,2) = w2016.teffZeuto1000.ave(6) - w2016.teffZeuto1000.stdev(6); % NA
metricsData(iB,iTeff,iW2016,1,1) = w2016.teffZeuto1000.ave(7);                                % STA, satellite export
metricsData(iB,iTeff,iW2016,2,1) = w2016.teffZeuto1000.ave(7) + w2016.teffZeuto1000.stdev(7); % STA
metricsData(iB,iTeff,iW2016,3,1) = w2016.teffZeuto1000.ave(7) - w2016.teffZeuto1000.stdev(7); % STA
metricsData(iB,iTeff,iW2016,1,2) = w2016.teffZeuto1000.ave(8);                                % STA, CMIP export
metricsData(iB,iTeff,iW2016,2,2) = w2016.teffZeuto1000.ave(8) + w2016.teffZeuto1000.stdev(8); % STA
metricsData(iB,iTeff,iW2016,3,2) = w2016.teffZeuto1000.ave(8) - w2016.teffZeuto1000.stdev(8); % STA
metricsData(iHo,iTeff,iW2016,1,1) = w2016.teffZeuto1000.ave(9);                                  % STP, satellite export
metricsData(iHo,iTeff,iW2016,2,1) = w2016.teffZeuto1000.ave(9) + w2016.teffZeuto1000.stdev(9);   % STP
metricsData(iHo,iTeff,iW2016,3,1) = w2016.teffZeuto1000.ave(9) - w2016.teffZeuto1000.stdev(9);   % STP
metricsData(iHo,iTeff,iW2016,1,2) = w2016.teffZeuto1000.ave(10);                                 % STP, CMIP export
metricsData(iHo,iTeff,iW2016,2,2) = w2016.teffZeuto1000.ave(10) + w2016.teffZeuto1000.stdev(10); % STP
metricsData(iHo,iTeff,iW2016,3,2) = w2016.teffZeuto1000.ave(10) - w2016.teffZeuto1000.stdev(10); % STP

ci_martinb = NaN(2,length(w2016.teffZeuto1000.ave));
midval_martinb = NaN(1,length(w2016.teffZeuto1000.ave));
ci_zstar = NaN(2,length(w2016.teffZeuto1000.ave));
midval_zstar = NaN(1,length(w2016.teffZeuto1000.ave));

for i = 1:length(w2016.teffZeuto1000.ave)
    % Error propagation using MC, assume gaussian distrib of errors
    A = generateMCparameters('gaussian',[w2016.teffZeuto1000.ave(i),w2016.teffZeuto1000.stdev(i)],'plot',false);
    A(A<0) = min(A(A>0));
    [midval_martinb(i),ci_martinb(:,i),funvals_martinb] = propagateErrorWithMC(funcMartinb_teffDep,A,'plot',false);  
    [midval_zstar(i),ci_zstar(:,i),funvals_zstar] = propagateErrorWithMC(funcZstar_teffDep,A,'plot',false);
end

metricsData(iE,iMartinb,iW2016,1,1) = midval_martinb(1);
metricsData(iE,iMartinb,iW2016,2,1) = ci_martinb(2,1);
metricsData(iE,iMartinb,iW2016,3,1) = ci_martinb(1,1);
metricsData(iE,iMartinb,iW2016,1,2) = midval_martinb(2);
metricsData(iE,iMartinb,iW2016,2,2) = ci_martinb(2,2);
metricsData(iE,iMartinb,iW2016,3,2) = ci_martinb(1,2);
metricsData(iO,iMartinb,iW2016,1,1) = midval_martinb(3);
metricsData(iO,iMartinb,iW2016,2,1) = ci_martinb(2,3);
metricsData(iO,iMartinb,iW2016,3,1) = ci_martinb(1,3);
metricsData(iO,iMartinb,iW2016,1,2) = midval_martinb(4);
metricsData(iO,iMartinb,iW2016,2,2) = ci_martinb(2,4);
metricsData(iO,iMartinb,iW2016,3,2) = ci_martinb(1,4);
metricsData(iP,iMartinb,iW2016,1,1) = midval_martinb(5);
metricsData(iP,iMartinb,iW2016,2,1) = ci_martinb(2,5);
metricsData(iP,iMartinb,iW2016,3,1) = ci_martinb(1,5);
metricsData(iP,iMartinb,iW2016,1,2) = midval_martinb(6);
metricsData(iP,iMartinb,iW2016,2,2) = ci_martinb(2,6);
metricsData(iP,iMartinb,iW2016,3,2) = ci_martinb(1,6);
metricsData(iB,iMartinb,iW2016,1,1) = midval_martinb(7);
metricsData(iB,iMartinb,iW2016,2,1) = ci_martinb(2,7);
metricsData(iB,iMartinb,iW2016,3,1) = ci_martinb(1,7);
metricsData(iB,iMartinb,iW2016,1,2) = midval_martinb(8);
metricsData(iB,iMartinb,iW2016,2,2) = ci_martinb(2,8);
metricsData(iB,iMartinb,iW2016,3,2) = ci_martinb(1,8);
metricsData(iHo,iMartinb,iW2016,1,1) = midval_martinb(9);
metricsData(iHo,iMartinb,iW2016,2,1) = ci_martinb(2,9);
metricsData(iHo,iMartinb,iW2016,3,1) = ci_martinb(1,9);
metricsData(iHo,iMartinb,iW2016,1,2) = midval_martinb(10);
metricsData(iHo,iMartinb,iW2016,2,2) = ci_martinb(2,10);
metricsData(iHo,iMartinb,iW2016,3,2) = ci_martinb(1,10);

metricsData(iE,iZstar,iW2016,1,1) = midval_zstar(1);
metricsData(iE,iZstar,iW2016,2,1) = ci_zstar(2,1);
metricsData(iE,iZstar,iW2016,3,1) = ci_zstar(1,1);
metricsData(iE,iZstar,iW2016,1,2) = midval_zstar(2);
metricsData(iE,iZstar,iW2016,2,2) = ci_zstar(2,2);
metricsData(iE,iZstar,iW2016,3,2) = ci_zstar(1,2);
metricsData(iO,iZstar,iW2016,1,1) = midval_zstar(3);
metricsData(iO,iZstar,iW2016,2,1) = ci_zstar(2,3);
metricsData(iO,iZstar,iW2016,3,1) = ci_zstar(1,3);
metricsData(iO,iZstar,iW2016,1,2) = midval_zstar(4);
metricsData(iO,iZstar,iW2016,2,2) = ci_zstar(2,4);
metricsData(iO,iZstar,iW2016,3,2) = ci_zstar(1,4);
metricsData(iP,iZstar,iW2016,1,1) = midval_zstar(5);
metricsData(iP,iZstar,iW2016,2,1) = ci_zstar(2,5);
metricsData(iP,iZstar,iW2016,3,1) = ci_zstar(1,5);
metricsData(iP,iZstar,iW2016,1,2) = midval_zstar(6);
metricsData(iP,iZstar,iW2016,2,2) = ci_zstar(2,6);
metricsData(iP,iZstar,iW2016,3,2) = ci_zstar(1,6);
metricsData(iB,iZstar,iW2016,1,1) = midval_zstar(7);
metricsData(iB,iZstar,iW2016,2,1) = ci_zstar(2,7);
metricsData(iB,iZstar,iW2016,3,1) = ci_zstar(1,7);
metricsData(iB,iZstar,iW2016,1,2) = midval_zstar(8);
metricsData(iB,iZstar,iW2016,2,2) = ci_zstar(2,8);
metricsData(iB,iZstar,iW2016,3,2) = ci_zstar(1,8);
metricsData(iHo,iZstar,iW2016,1,1) = midval_zstar(9);
metricsData(iHo,iZstar,iW2016,2,1) = ci_zstar(2,9);
metricsData(iHo,iZstar,iW2016,3,1) = ci_zstar(1,9);
metricsData(iHo,iZstar,iW2016,1,2) = midval_zstar(10);
metricsData(iHo,iZstar,iW2016,2,2) = ci_zstar(2,10);
metricsData(iHo,iZstar,iW2016,3,2) = ci_zstar(1,10);

fprintf('...done.\n')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 7 - SAVE THE DATA
% -------------------------------------------------------------------------

save(fullfile('.','data','processed',filenameOutputMetricsData),'metricsData',...
    'nPublications','nMetrics')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 8 - WRITE OUT THE METRICS DATA 5D ARRAY INTO A .CSV FILE 
% (DATASET S2)
% -------------------------------------------------------------------------

% I need to flatten or reshape the array into a 2D matrix 

pubNames = {'Francois et al. (2002)',...
            'Buesseler & Boyd (2009)',...
            'Lam et al. (2011)',...
            'Guidi et al. (2015)',...
            'Mouw et al. (2016b)',...
            'Henson et al. (2012)',...                  
            'Marsay et al. (2015)',...
            'Weber et al. (2016)',...
            'UVP5 compilation (Kiko et al. 2022)',...
            'Trap & radionuclide compilation (this study)'}';

% Quadruplicate each entry
quadruplicateRefNames = repelem(pubNames, 4); % up to 4 values provided by a publication

% Initialise the output array
outputArray = NaN(length(quadruplicateRefNames),(3*NUM_LOCS*3));

%                                     Martin's b                                   z*    Teff  
%        ----------------------------------------------------------------------- ------ ------
%         HOT/ALOHA    BATS/OFP      EqPac      PAP-SO       OSP      HAUSGARTEN   ...    ...
%        ----------- ----------- ----------- ----------- ----------- ----------- 
%        -std m +std -std m +std -std m +std -std m +std -std m +std -std m +std   ...    ...
%        --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
% Ref     .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .    ...    ...   
% Ref
% ...

funzstar = @(x) sprintf('%3.0f', x);
funrest = @(x) sprintf('%.3f', x);

iRow = 0; % index for the output array
iColumn = 0; % index for the output array
for iRef = 1:nPublications
    for iRep = 1:4
        iRow = iRow + 1;
        iColumn = 0;
        for iMetric = [iMartinb,iZstar,iTeff]
            for iLoc = [iHo,iB,iE,iP,iO,iHa]
                for iStatistic = [3,1,2] % low std, mean, upp std
                    iColumn = iColumn + 1;
                    extract = metricsData(iLoc,iMetric,iRef,iStatistic,iRep);
                    % Cap the number of decimal digits accordingly
                    if (iMetric == iZstar && ~isnan(extract))
                        C = num2cell(extract);
                        F = cellfun(funzstar, C, 'UniformOutput',0);
                        extract = str2double(F);
                    elseif (iMetric ~= iZstar && ~isnan(extract))
                        C = num2cell(extract);
                        F = cellfun(funrest, C, 'UniformOutput',0);
                        extract = str2double(F);
                    end
                    outputArray(iRow,iColumn) = extract;
                end
            end % iLoc
        end % iMetric
    end % iRep
end % iRef

% Delete rows where all values are NaN
nanRows = all(isnan(outputArray), 2);
outputArray(nanRows,:) = [];
quadruplicateRefNames(nanRows) = [];

% Concatenate row names and numeric array horizontally...
horzConcatenatedArray = [quadruplicateRefNames,num2cell(outputArray)];

% ... and now vertically with the headers
headerRowFirst = [{''},repmat({'Martin_b'},1,(NUM_LOCS*3)),...
    repmat({'z_star'},1,(NUM_LOCS*3)),repmat({'Teff_100_to_1000m'},1,(NUM_LOCS*3))];

headerRowSecond = [{''},repmat([repmat({'HOT/ALOHA'},1,3),repmat({'BATS/OFP'},1,3),repmat({'EqPac'},1,3),...
    repmat({'PAP-SO'},1,3),repmat({'OSP'},1,3),repmat({'HAUSGARTEN'},1,3)],1,3)];

headerRowThird = ['Reference',repmat({'low_std','mean','upp_std'},1,NUM_LOCS*3)];

outputTable = cell2table([headerRowFirst;...
                          headerRowSecond;...
                          headerRowThird;...
                          horzConcatenatedArray]);

% Write table to CSV without column names
writetable(outputTable,fullfile('.','data','processed',filenameOutputCsvTable),...
    'WriteVariableNames',false,'Delimiter',',');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 9 - ANOVA STATISTICAL TESTS
% -------------------------------------------------------------------------

% ANOVA (Analysis of Variance) and the Kruskal-Wallis test are statistical 
% tests used to compare groups, but they are applicable under different
% conditions: ANOVA is used when the data points are normally distributed
% and variances are homogeneous, and Kruskal-Wallis when that is not the 
% case.
%
% We want to answer these two questions (run two separate tests):
% (1) Do locations have metric values that are statistically different?
% (2) Do publications produce metric values that are statistically 
% different? Said otherwise, can we see a "publication" footprint in the
% metrics?

% .........................................................................

load(fullfile('.','data','processed',filenameOutputMetricsData),'metricsData',...
    'nPublications','nMetrics')

%% ANOVA test to see the effect of geographical location on the different metrics

fprintf('\nANOVA test to see the effect of geographical location on the different BCP metrics\n')

% Put all samples from each reference together
localMetricSamples = NaN(3,NUM_LOCS,nPublications*4); % 4 for the number of repetitions
for iMetric = [iMartinb,iZstar,iTeff] % response variable
    for iLoc = [iHo,iB,iE,iP,iO,iHa] % treatment (independent variable)
        samples = squeeze(metricsData(iLoc,iMetric,:,1,:)); % 1 for the mean value
        localMetricSamples(iMetric,iLoc,:) = reshape(samples,1,[]);
    end
end

pvalueTableShapiroTest = zeros(3,NUM_LOCS); % 3 metrics x nLocs
pValueBartlettTest = zeros(3,1);
pValue = zeros(3,1);

i = 0;
for iMetric = [iMartinb,iZstar,iTeff] % response variable
    
    fprintf('This is metric %0.f\n',iMetric)
    i = i + 1;
    
    % Factor "location" with 6 levels (k=6)
    % N=10 replicates/level (nReferences)
    loc1 = squeeze(localMetricSamples(iMetric,1,:)); % HOT/ALOHA
    loc2 = squeeze(localMetricSamples(iMetric,2,:)); % BATS/OFP
    loc3 = squeeze(localMetricSamples(iMetric,3,:)); % EqPac 
    loc4 = squeeze(localMetricSamples(iMetric,4,:)); % PAP-SO
    loc5 = squeeze(localMetricSamples(iMetric,5,:)); % OSP
    loc6 = squeeze(localMetricSamples(iMetric,6,:)); % HAUSGARTEN

    data = [loc1, loc2, loc3, loc4, loc5, loc6]; % each column is a group (level, treatment, location)

    % Testing if data are normally distributed -Shapiro-Wilk test
    for iGroup = 1:NUM_LOCS
        [H, pvalueTableShapiroTest(i,iGroup), W] = swtest(data(:,iGroup)); %adtest(data(:,iCol));
    end

    % Testing if variances are homogeneous -Bartlett test
    pValueBartlettTest(i) = vartestn(data,'TestType','Bartlett','Display','off');
    
    % Decision tree
    if all(pvalueTableShapiroTest(i,:) > 0.05, 'all') && pValueBartlettTest(i) > 0.05
        % Use ANOVA
        disp("Use ANOVA")
        [pValue(i),tbl,stats] = anova1(data);
    else
        % Use Kruskal-Wallis
        disp("Use Kruskal-Wallis")
        [pValue(i),tbl,stats] = kruskalwallis(data);
    end
    
    % If the test rejects the null hypothesis that all group means are 
    % equal, we can use the multiple comparisons to determine which group 
    % means are different from others. For that, we will use multcompare.
    if (pValue(i) < 0.05) % rejects H0
        disp('We reject H0 that all group means are equal')
        [c,m,h,nms] = multcompare(stats); % , 'Display', 'off'
    else
        disp('We accept H0 that all group means are equal')
    end
        
end

% .........................................................................

%% ANOVA test to see the effect of publication on the different BCP metrics

fprintf('\nANOVA test to see the effect of publication on the different BCP metrics\n')

% Put all samples from each reference together
refMetricSamples = NaN(3,nPublications,NUM_LOCS*4); % 4 for the number of repetitions
for iMetric = [iMartinb,iZstar,iTeff] % response variable
    for iReference = 1:nPublications % treatment (independent variable)
        samples = squeeze(metricsData(:,iMetric,iReference,1,:)); % 1 for the mean value
        refMetricSamples(iMetric,iReference,:) = reshape(samples,1,[]);
    end
end

pvalueTableShapiroTest = zeros(3,nPublications); % 3 metrics x nReferences
pValueBartlettTest = zeros(3,1);
pValue = zeros(3,1);

i = 0;
for iMetric = [iMartinb,iZstar,iTeff] % response variable
    
    fprintf('This is metric %0.f\n',iMetric)
    i = i + 1;
    
    % Factor "reference" with 10 levels (k=10)
    % N=6*4 replicates/level (nLocs*nRepetitions)
    ref1 = squeeze(refMetricSamples(iMetric,1,:)); 
    ref2 = squeeze(refMetricSamples(iMetric,2,:)); 
    ref3 = squeeze(refMetricSamples(iMetric,3,:)); 
    ref4 = squeeze(refMetricSamples(iMetric,4,:));
    ref5 = squeeze(refMetricSamples(iMetric,5,:));
    ref6 = squeeze(refMetricSamples(iMetric,6,:));
    ref7 = squeeze(refMetricSamples(iMetric,7,:));
    ref8 = squeeze(refMetricSamples(iMetric,8,:));
    ref9 = squeeze(refMetricSamples(iMetric,9,:));
    ref10 = squeeze(refMetricSamples(iMetric,10,:));

    data = [ref1, ref2, ref3, ref4, ref5,...
            ref6, ref7, ref8, ref9, ref10]; % each column is a group (level, treatment, reference)

    % Testing if data are normally distributed -Shapiro-Wilk test
    for iGroup = 1:nPublications
        [H, pvalueTableShapiroTest(i,iGroup), W] = swtest(data(:,iGroup)); %adtest(data(:,iCol));
    end

    % Testing if variances are homogeneous -Bartlett test
    pValueBartlettTest(i) = vartestn(data,'TestType','Bartlett','Display','off');
    
    % Decision tree
    if all(pvalueTableShapiroTest(i,:) > 0.05, 'all') && pValueBartlettTest(i) > 0.05
        % Use ANOVA
        disp("Use ANOVA")
        [pValue(i),tbl,stats] = anova1(data);
    else
        % Use Kruskal-Wallis
        disp("Use Kruskal-Wallis")
        [pValue(i),tbl,stats] = kruskalwallis(data);
    end
    
    % If the test rejects the null hypothesis that all group means are 
    % equal, we can use the multiple comparisons to determine which group 
    % means are different from others. For that, we will use multcompare
    % (Tukey-Kramer post-hoc test).
    if (pValue(i) < 0.05) % rejects H0
        disp('We reject H0 that all group means are equal')
        [c,m,h,nms] = multcompare(stats); % , 'Display', 'off'
    else
        disp('We accept H0 that all group means are equal')
    end
        
end
