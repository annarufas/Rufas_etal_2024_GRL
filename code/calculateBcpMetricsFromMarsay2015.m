function [marsay2015] = calculateBcpMetricsFromMarsay2015(...
    isGlobalCalculation,listLocalLats,listLocalLons)

% CALCULATEBCPMETRICSFROMMARSAY2015 Calculates three BCP metrics: Martin's
% b and z* using the statistical fit of Marsay et al. (2015), as well Teff 
% as a function of b. The three metrics are calculated using the annual 
% means of the variables they have dependencies on.
%
%   INPUT: 
%       isGlobalCalculation - choice (1 or 0)
%       listLocalLats       - list of local latitudes
%       listLocalLons       - list of local longitudes
%
%   OUTPUT:
%       isGlobalCalculation = 1, useful for visualising the global output 
%           of the Marsay et al. (2015) algorithm. The output is a global 
%           array of annual mean values of the three BCP metrics, with no 
%           standard deviations.
%       isGlobalCalculation = 0 produces local annual mean values and 
%           standard deviations of the three BCP metrics.
% 
%   This script uses these external functions: 
%       generateMCparameters.m - from FileExchange
%       propagateErrorWithMC.m - from FileExchange
% 
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 5 Nov 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Define parameters of the fitted b-T and z-T curves

% Mean values for 1=slope b-T, 2=intercept b-T, 3=slope z-T, 4=intercept z-T
fitParam.mid = [0.062, 0.303, 19.0, 483]; 

%% Calculate the std of the fitted parameters in Marsay (local function)

% This function also recalculates the mean values shown above and since
% they ar, almost identical to the ones provided by Marsay et al. (2015), but
% not exactly the same, I have updated them.
[fmartinb,fzstar] = calculateMarsay2015fitParams(fitParam.mid); 

fitParam.ave   = [fmartinb.slope.ave,   fmartinb.intercept.ave,   fzstar.slope.ave,    fzstar.intercept.ave]; % add updated mean values
fitParam.stdev = [fmartinb.slope.stdev, fmartinb.intercept.stdev, fzstar.slope.stdev,  fzstar.intercept.stdev];

%% Define Marsay algorithm (Eqs. in Fig. 2a and 2b in Marsay et al. (2015)) 

funcMartinb = @(x) (x(1)*x(2)) + x(3); % b = 0.062.*tempUp500 + 0.303 
funcZstar   = @(x) x(1) - (x(2)*x(3)); % z* = 483 - 19.*tempUp500 

%% Define supplementary function for Teff 100 to 1000 m

z0 = 100;
z1 = 1000;
funcTeff = @(x) (z1/z0).^(-x); % 'x' is Martin b

%% Load data inputs for Marsay's algorithm

load(fullfile('.','data','interim','temp_annual_woa23.mat'),...
    'temp_annual','temp_annual_std','woa_annual_lat','woa_annual_lon','woa_annual_depth_temp')

% Extract temperature in the first 500 m of the water column and calculate
% the median and propagate errors
iDepth500 = find(woa_annual_depth_temp == 500);

tempAnnualUpp500_vals = NaN(size(temp_annual,1),size(temp_annual,2));
tempAnnualUpp500_err = NaN(size(tempAnnualUpp500_vals));

for iLat = 1:length(woa_annual_lat)
    for iLon = 1:length(woa_annual_lon)

        ta = squeeze(temp_annual(iLat,iLon,(1:iDepth500)));
        td = squeeze(temp_annual_std(iLat,iLon,(1:iDepth500)));
        idxNonNan = ~isnan(ta) & ~isnan(td);
        ta = ta(idxNonNan);
        td = td(idxNonNan);

        % If there are valid data points, calculate median and propagate error
        if (any(ta))
            % Calculate the median of the temperature values
            tempAnnualUpp500_vals(iLat,iLon) = median(ta);

            % Error propagation using worstcase with the function handle and the parameters
            [~,~,~,f_MID,f_UB,~,~] = worstcase(@(ta) median(ta),ta,td);
            tempAnnualUpp500_err(iLat,iLon) = f_UB - f_MID;
        end
    end
end

%% Calculations

if isGlobalCalculation

    % Initialise output arrays
    [nr_temp,nc_temp] = size(tempAnnualUpp500_median);
    martinb = NaN(nr_temp,nc_temp);
    zstar = NaN(size(martinb));
    
    % Calculate Martin's b and z*
    for iRow = 1:nr_temp
        for iCol = 1:nc_temp
            temp = tempAnnualUpp500_vals(iRow,iCol);
            if (~isnan(temp))
                paramArrayMartinb = [fitParam.ave(1);temp;fitParam.ave(2)];
                paramArrayZstar = [fitParam.ave(4);fitParam.ave(3);temp];
                martinb(iRow,iCol) = funcMartinb(paramArrayMartinb);
                zstar(iRow,iCol) = funcZstar(paramArrayZstar);
            end
        end
    end
    
    % Calculate Teff 100 to 1000 m as a function of b
    teff1000 = arrayfun(funcTeff,martinb);

    % Save arrays
    marsay2015.martinb  = martinb;
    marsay2015.zstar    = zstar;
    marsay2015.teff1000 = teff1000;
    
    save(fullfile('.','data','processed','bcpmetrics_marsay2015_global.mat'),'marsay2015*')

%     % Visual checks
%     figure(); pcolor(flipud(rot90(tempAnnualUpp500(:,:)))); caxis([-2 25]); 
%     colorbar; colormap(jet); shading interp; title('WOA23 median temperature upp. 500 m')
%     figure(); pcolor(flipud(rot90(martinb(:,:)))); caxis([0.26 2]); 
%     colorbar; colormap(jet); shading interp; title('Marsay et al. (2015), b')
%     figure(); pcolor(flipud(rot90(zstar(:,:)))); caxis([0 600]); 
%     colorbar; colormap(jet); shading interp; title('Marsay et al. (2015), z* (m)')
%     figure(); pcolor(flipud(rot90(teff1000(:,:)))); caxis([0 0.4]); 
%     colorbar; colormap(jet); shading interp; title('Marsay et al. (2015), Teff 100 to 1000 m')
  
else % local data
    
    qLons = listLocalLons;   
    qLats = listLocalLats;
    nLocs = length(listLocalLats);
    
    % Initialise output arrays
    trueMartinb  = NaN(nLocs,1);
    trueZstar    = NaN(size(trueMartinb));
    estimMartinb = NaN(size(trueMartinb));
    estimZstar   = NaN(size(trueMartinb));
    
    % Prepare interpolant to extract local temperature for our locations of interest
    [X, Y] = ndgrid(woa_annual_lat, woa_annual_lon); 
    Fv = griddedInterpolant(X, Y, tempAnnualUpp500_vals, 'linear');
    Fe = griddedInterpolant(X, Y, tempAnnualUpp500_err, 'linear');

    for iLoc = 1:nLocs

        [qX, qY] = ndgrid(qLats(iLoc),qLons(iLoc)); % query points for interpolation
        qTempAnnualMean = Fv(qX, qY);
        qTempAnnualStd = Fe(qX, qY);
        
        % Calculate Martin's b and z*
        paramArray = [fitParam.ave(1);qTempAnnualMean;fitParam.ave(2)];
        trueMartinb(iLoc) = funcMartinb(paramArray); 
        paramArray = [fitParam.ave(4);fitParam.ave(3);qTempAnnualMean];
        trueZstar(iLoc) = funcZstar(paramArray);
        
        % Error propagation using MC
        A = generateMCparameters('gaussian',[fitParam.ave(1),fitParam.stdev(1)],'plot',false);
        B = generateMCparameters('gaussian',[qTempAnnualMean,qTempAnnualStd],'plot',false);
        C = generateMCparameters('gaussian',[fitParam.ave(2),fitParam.stdev(2)],'plot',false);
        [midval_martinb,ci_martinb,funvals_martinb] = propagateErrorWithMC(funcMartinb,[A;B;C],'plot',false);  
        
        A = generateMCparameters('gaussian',[fitParam.ave(4),fitParam.stdev(4)],'plot',false);
        B = generateMCparameters('gaussian',[fitParam.ave(3),fitParam.stdev(3)],'plot',false);
        C = generateMCparameters('gaussian',[qTempAnnualMean,qTempAnnualStd],'plot',false);
        [midval_zstar,ci_zstar,funvals_zstar] = propagateErrorWithMC(funcZstar,[A;B;C],'plot',false);
        midval_zstar(midval_zstar < 0) = 0;
        ci_zstar(ci_zstar < 0) = 0;

        estimMartinb(iLoc) = midval_martinb;
        estimZstar(iLoc)   = midval_zstar;

        % Calculate Teff 100 to 1000 m as a function of b 
        if (~isnan(midval_martinb))
            
            % The following three methods produce similar results. I am
            % going to use Method C as it requires less steps
            
%             % Method A
%             A = generateMCparameters('bootstrapDistribution',funvals_martinb);
%             [midval_teff1000_1,ci_teff1000_1,funvals_teff1000_1] = propagateErrorWithMC(funcTeff,A,'plot',false);
% 
%             % Method B - caveat: the error distribution is not gaussian
%             B = generateMCparameters('gaussian',[midval_martinb,(ci_martinb(2)-ci_martinb(1))/2]);
%             [midval_teff1000_2,ci_teff1000_2,funvals_teff1000_2] = propagateErrorWithMC(funcTeff,B,'plot',false);
%             
%             % Method C
%             [midval_teff1000_3,ci_teff1000_3,funvals_teff1000_3] = propagateErrorWithMC(funcTeff,funvals_martinb,'plot',false);

            [midval_teff1000,ci_teff1000,funvals_teff1000] =...
                propagateErrorWithMC(funcTeff,funvals_martinb,'plot',false);
            
        end

        if (nLocs < 10)
            fprintf('\nLocation %d',iLoc)
            fprintf('\nThe local T is %4.2f and the std is %4.2f.',qTempAnnualMean,qTempAnnualStd)
            fprintf('\nThe true b is %4.3f and the estimated b is %4.3f.',trueMartinb(iLoc),estimMartinb(iLoc))
            fprintf('\nThe bounds estimated for b are %4.3f to %4.3f.',ci_martinb(1),ci_martinb(2))
            fprintf('\nThe true z* is %3.1f and the estimated z* is %3.1f.',trueZstar(iLoc),estimZstar(iLoc))
            fprintf('\nThe bounds estimated for z* are %3.1f to %3.1f.',ci_zstar(1),ci_zstar(2))
%             fprintf('\nThe estimated Teff 100 to 1000 m Method A is %4.3f, with bounds %4.3f to %4.3f',midval_teff1000_1,ci_teff1000_1(1),ci_teff1000_1(2))
%             fprintf('\nThe estimated Teff 100 to 1000 m Method B is %4.3f, with bounds %4.3f to %4.3f',midval_teff1000_2,ci_teff1000_2(1),ci_teff1000_2(2))
%             fprintf('\nThe estimated Teff 100 to 1000 m Method C is %4.3f, with bounds %4.3f to %4.3f\n',midval_teff1000_3,ci_teff1000_3(1),ci_teff1000_3(2))
            fprintf('\nThe estimated Teff is %4.3f, with bounds %4.3f to %4.3f.\n',midval_teff1000,ci_teff1000(1),ci_teff1000(2))
        end
        
        % Save arrays
        marsay2015.martinb.ave(iLoc)       = midval_martinb;
        marsay2015.martinb.stdevupp(iLoc)  = ci_martinb(2);
        marsay2015.martinb.stdevlow(iLoc)  = ci_martinb(1);
        marsay2015.martinb.max(iLoc)       = max(funvals_martinb);
        marsay2015.martinb.min(iLoc)       = min(funvals_martinb);
        
        marsay2015.zstar.ave(iLoc)         = midval_zstar;
        marsay2015.zstar.stdevupp(iLoc)    = ci_zstar(2);
        marsay2015.zstar.stdevlow(iLoc)    = ci_zstar(1);
        marsay2015.zstar.max(iLoc)         = max(funvals_zstar);
        marsay2015.zstar.min(iLoc)         = min(funvals_zstar);
        
        marsay2015.teff1000.ave(iLoc)      = midval_teff1000;
        marsay2015.teff1000.stdevupp(iLoc) = ci_teff1000(2);
        marsay2015.teff1000.stdevlow(iLoc) = ci_teff1000(1);
        marsay2015.teff1000.max(iLoc)      = max(funvals_teff1000);
        marsay2015.teff1000.min(iLoc)      = min(funvals_teff1000);

    end % iLoc

    save(fullfile('.','data','processed','bcpmetrics_marsay2015_local.mat'),'marsay2015*')

end % isGriddedData

%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

% *************************************************************************

function [fmartinb,fzstar] = calculateMarsay2015fitParams(marsayParams)

% CALCULATEMARSAY2015FITPARAMS Recalculates the fit parameters of the 
% equations z*-T and b-T using the data points that Marsay et al. (2015) 
% report in their Suppl. This is to check that the equation that Marsay & Co
% typed in their Fig. 2b does not contain a typo (after I realised that z* 
% is too low in the subtropical and tropical oceans (<50 m)). This check 
% also provides the opportunity to estimate the standard deviation around 
% the mean of the fit parameters. Whereas z* and b values are provided in 
% the Suppl., temperature values are not, so we've digitised Marsay's 
% Fig. 2b to extract temperature. We've also checked that those temperature 
% values make sense by extracting them from the WOA18 monthly climatology 
% and calculating the median temperature for the first 500 m of the water 
% column for the months in which POC flux was sampled to calculate z*. 

% Calculate the mean lat and lon for the NASG region
ts1.nasg.lat = [ 26.27,  26.60,  26.82,  26.27,  25.66,  26.27,  26.82,  26.60,  25.66,  26.20];
ts1.nasg.lon = [-31.08, -31.33, -31.59, -31.08, -30.32, -31.08, -31.59, -31.33, -30.32, -31.32];

% The following data are from Table S3
%
%              PAP, Irminger, Iceland,               NASG, K2-1, K2-2, ALOHA-1, ALOHA-2
ts3.lat     = [ 48.92,  60.84,  62.15, mean(ts1.nasg.lat),   47,   47,  22.5,  22.5];
ts3.lon     = [-16.39, -31.60, -24.39, mean(ts1.nasg.lon),  160,  160,  -158,  -158];
ts3.zstar   = [   259,    234,    296,                135,  442,  523,   223,   177];
ts3.martinb = [  0.70,   0.88,   0.69,               1.59, 0.57, 0.49,  1.25,  1.36];
ts3.temp    = [  11.6,    6.3,    7.9,               17.7,  3.2,  3.2,  16.6,  15.8]; % WebPlotDigitizer
ts3.month   = [     8,      7,      8,                  8,    7,    8,     6,     7];

x = ts3.temp'; % temperature vector
yb = ts3.martinb';
yz = ts3.zstar';
n = length(x); % no. samples

% .........................................................................

% Fit a regression line between b and temperature
fitMartinb = fittype('mb.*x + nb','independent','x','dependent','yb','coefficients',{'mb','nb'});
[fb,~] = fit(x,yb,fitMartinb,'StartPoint',[0.05,0.1]);

% The fitted parameters are
fmartinb.slope.ave = fb.mb;
fmartinb.intercept.ave = fb.nb;

% Extract the 95% CI bounds of the coefficients m and n
coeffBounds = confint(fb); 
fmartinb.slope.CI95 = coeffBounds(:,1);
fmartinb.intercept.CI95 = coeffBounds(:,2);

% Calculate std from 95% CI
fmartinb.slope.stdev = ((max(fmartinb.slope.CI95)-min(fmartinb.slope.CI95))*sqrt(n))/3.92; 
fmartinb.intercept.stdev = ((max(fmartinb.intercept.CI95)-min(fmartinb.intercept.CI95))*sqrt(n))/3.92;

fprintf('\nb-T eqn.')
fprintf('\nThe slope calculated by Marsay is %4.3f,',marsayParams(1))
fprintf('\nand the slope calculated here is %4.3f,',fmartinb.slope.ave)
fprintf('\nwith 95CI boundaries %4.3f to %4.3f and std %4.3f.',max(fmartinb.slope.CI95),min(fmartinb.intercept.CI95),fmartinb.slope.stdev)
fprintf('\nThe intercept calculated by Marsay is %4.3f,',marsayParams(2))
fprintf('\nand the intercept calculated here is %4.3f,',fmartinb.intercept.ave)
fprintf('\nwith 95CI boundaries %4.3f to %4.3f and std %4.3f.\n',max(fmartinb.intercept.CI95),min(fmartinb.intercept.CI95),fmartinb.intercept.stdev)

%figure(); plot(fb,x,yb); xlim([0 20]); ylim([0 2])
%xlabel('Temperature (ºC)'); ylabel('Martin b'); title('Check Marsay algorithm');

% .........................................................................

% Fit a regression line between z* and temperature
fitZstar = fittype('nz - mz.*x','independent','x','dependent','yz','coefficients',{'mz','nz'});
[fz,~] = fit(x,yz,fitZstar,'StartPoint',[400,10]);

% The fitted parameters are
fzstar.slope.ave = fz.mz;
fzstar.intercept.ave = fz.nz;

% Extract the 95% CI bounds of the coefficients m and n
coeffBounds = confint(fz); 
fzstar.slope.CI95 = coeffBounds(:,1);
fzstar.intercept.CI95 = coeffBounds(:,2);

% Calculate std from 95% CI
fzstar.slope.stdev = ((max(fzstar.slope.CI95)-min(fzstar.slope.CI95))*sqrt(n))/3.92; 
fzstar.intercept.stdev = ((max(fzstar.intercept.CI95)-min(fzstar.intercept.CI95))*sqrt(n))/3.92;

fprintf('\nz-T eqn.')
fprintf('\nThe slope calculated by Marsay is %3.0f,',marsayParams(3))
fprintf('\nand the slope calculated here is %3.0f,',fzstar.slope.ave)
fprintf('\nwith 95CI boundaries %3.0f to %3.0f and std %3.0f.',max(fzstar.slope.CI95),min(fzstar.slope.CI95),fzstar.slope.stdev)
fprintf('\nThe intercept calculated by Marsay is %3.0f,',marsayParams(4))
fprintf('\nand the intercept calculated here is %3.0f,',fzstar.intercept.ave)
fprintf('\nwith 95CI boundaries %3.0f to %3.0f and std %3.0f.\n',max(fzstar.intercept.CI95),min(fzstar.intercept.CI95),fzstar.intercept.stdev)

%figure(); plot(fz,x,yz); xlim([0 20]); ylim([0 800])
%xlabel('Temperature (ºC)'); ylabel('z^{*}'); title('Check Marsay algorithm'); 

end % calculateMarsay2015fitParams

% *************************************************************************

end