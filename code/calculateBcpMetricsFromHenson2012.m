function [henson2012] = calculateBcpMetricsFromHenson2012(...
    isGlobalCalculation,listLocalLats,listLocalLons)

% CALCULATEBCPMETRICSFROMHENSON2012 Calculates three BCP metrics: b and 
% z* and Teff as a function of b. The three metrics are calculated using 
% the annual means of the variables they have dependencies on.
%
%   INPUT: 
%       isGlobalCalculation - choice (1 or 0)
%       listLocalLats       - list of local latitudes
%       listLocalLons       - list of local longitudes
%
%   OUTPUT:
%       isGlobalCalculation = 1, useful for visualising the global output 
%           of the Henson et al. (2012) algorithm. The output is a global 
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
%   Version 1.0 - Completed 4 Nov 2024 
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Define parameters of the fitted PEeff-T curve

% PEeff = 0.23*exp(-0.08*SST), where the std is provided in the caption of 
% Fig. 1 in Henson et al. (2011)

fitParam.ave   = [0.23, 0.08];
fitParam.stdev = [0.04, 0.01]; 
fitParam.stdevlow = [fitParam.ave(1)-fitParam.stdev(1), fitParam.ave(2)-fitParam.stdev(2)];
fitParam.stdevupp = [fitParam.ave(1)+fitParam.stdev(1), fitParam.ave(2)+fitParam.stdev(2)];

%% Define supplementary functions for Teff 100 to 1000 m and z*

% Notice that Henson et al. (2012) provide Teff from 100 to 200 m, so Teff 
% from 100 to 1000 m needs to be calculated, and we use Martin's b for
% that.

z0 = 100;
z1 = 1000;
funcTeff  = @(x) (z1/z0).^(-x); % 'x' is Martin b
funcZstar = @(x) -(z1-z0)./log((z1/z0).^(-x)); % 'x' is Martin b

%% Load NPP from Carr 2002 algorithm and SST (global arrays)

load(fullfile('.','data','interim','npp_carr2002_seawifs_pathfinder.mat'),...
    'npp_avg','npp_lat','npp_lon')
load(fullfile('.','data','interim','sst_pathfinder_v5.mat'),...
    'sst','sst_lat','sst_lon')

% Regrid SST array to NPP grid
[qX, qY, qT] = ndgrid(npp_lat, npp_lon, (1:12)'); % query points for interpolation    
[X, Y, T] = ndgrid(sst_lat, sst_lon, (1:12)'); % original array
F = griddedInterpolant(X, Y, T, sst, 'linear', 'none'); % interpolant
sst_regridded = F(qX, qY, qT); % regridded SST array
    
%% Calculations

% Apply Henson's algorithm to calculate b, PEeff and Teff 100 to 2000 m.
% No error propagation in the global calculations version (time-consuming).
     
if isGlobalCalculation

    % Prepare input data arrays
    qNppAnnualMean = mean(npp_avg,3,'omitnan'); % mg C m-2 d-1
    qNppAnnualStd = std(npp_avg,0,3,'omitnan');
    qSstAnnualMean = mean(sst_regridded,3,'omitnan'); % ºC

    % Initialise output arrays
    [nr_npp,nc_npp] = size(qNppAnnualMean);
    martinb  = NaN(nr_npp,nc_npp);
    teff2000 = NaN(size(martinb));
    peeff    = NaN(size(martinb));

    % Calculate PEeff, Teff 100 to 2000 m and b
    for iRow = 1:nr_npp
        for iCol = 1:nc_npp
            if (qNppAnnualMean(iRow,iCol) > 0 && ~isnan(qNppAnnualMean(iRow,iCol)))
                svi = qNppAnnualStd(iRow,iCol)./qNppAnnualMean(iRow,iCol); % seasonal variation index of NPP
                paramArray = [qSstAnnualMean(iRow,iCol);fitParam.ave(1);fitParam.ave(2);qNppAnnualMean(iRow,iCol);svi];
                peeff(iRow,iCol) = Henson2012peeff(paramArray);  
                teff2000(iRow,iCol) = Henson2012teff2000(paramArray);
                martinb(iRow,iCol) = Henson2012martinb(paramArray);
            end
        end
    end
    
    % Calculate z* and Teff 100 to 1000 m as a function of Martin's b
    zstar = arrayfun(funcZstar,martinb);
    teff1000 = arrayfun(funcTeff,martinb);

    % Save outputs
    henson2012.martinb  = martinb;
    henson2012.peeff    = peeff;
    henson2012.teff2000 = teff2000;
    henson2012.teff1000 = teff1000;
    henson2012.zstar    = zstar;

    save(fullfile('.','data','processed','bcpmetrics_henson2012_global.mat'),'henson2012*')
 
%     % Visual check
%     figure(); pcolor(flipud(rot90(martinb(:,:)))); caxis([0.2 1.2]);
%     colorbar; colormap(jet); shading interp; title('Henson 2012 Martin b')

else % local data
    
    qLons = listLocalLons;   
    qLats = listLocalLats;
    nLocs = length(listLocalLats);
    
    % Prepare interpolants to extract local NPP and SST for our locations of interest
    [X, Y, T] = ndgrid(npp_lat, npp_lon, (1:12)'); % same for SST
    Fnpp = griddedInterpolant(X, Y, T, npp_avg, 'linear');
    Fsst = griddedInterpolant(X, Y, T, sst_regridded, 'linear'); 

    % Initialise output arrays
    trueMartinb   = NaN(nLocs,1);
    trueTeff2000  = NaN(size(trueMartinb));
    truePeeff     = NaN(size(trueMartinb));
    estimMartinb  = NaN(size(trueMartinb));
    estimTeff2000 = NaN(size(trueMartinb));
    estimPeeff    = NaN(size(trueMartinb));

    for iLoc = 1:nLocs

        [qX, qY, qT] = ndgrid(qLats(iLoc), qLons(iLoc), (1:12)'); % query points for interpolation
        qNppMonthly = Fnpp(qX, qY, qT);
        qSstMonthly = Fsst(qX, qY, qT);

        % Prepare input data arrays
        qNppAnnualMean = mean(qNppMonthly,'omitnan'); % mg C m-2 d-1
        qNppAnnualStd = std(qNppMonthly,0,'omitnan');        
        qSstAnnualMean = mean(qSstMonthly,'omitnan'); % ºC
        qSstAnnualStd = std(qSstMonthly,0,'omitnan');
        
        % Calculate PEeff, Teff 100 to 2000 m and b
        if (qNppAnnualMean > 0 && ~isnan(qNppAnnualMean))
            svi = qNppAnnualStd/qNppAnnualMean; % seasonal variation index of NPP
            paramArray = [qSstAnnualMean;fitParam.ave(1);fitParam.ave(2);qNppAnnualMean;svi];
            truePeeff(iLoc) = Henson2012peeff(paramArray);  
            trueTeff2000(iLoc) = Henson2012teff2000(paramArray);
            trueMartinb(iLoc) = Henson2012martinb(paramArray);
        end
   
        % Error propagation using MC
        A = generateMCparameters('gaussian',[qSstAnnualMean,qSstAnnualStd],'plot',false);
        B = generateMCparameters('gaussian',[fitParam.ave(1),fitParam.stdev(1)],'plot',false);
        C = generateMCparameters('gaussian',[fitParam.ave(2),fitParam.stdev(2)],'plot',false);
        D = generateMCparameters('gaussian',[qNppAnnualMean,qNppAnnualStd],'plot',false);
        D(D<0) = min(D(D>0));
        E = repmat(svi,[1 length(A)]);
        paramMatrix = [A;B;C;D;E];
        [midval_peeff,ci_peeff,funvals_peeff] = propagateErrorWithMC(@Henson2012peeff,paramMatrix,'plot',false);  
        [midval_teff2000,ci_teff2000,funvals_teff2000] = propagateErrorWithMC(@Henson2012teff2000,paramMatrix,'plot',false);
        [midval_martinb,ci_martinb,funvals_martinb] = propagateErrorWithMC(@Henson2012martinb,paramMatrix,'plot',false);

        estimMartinb(iLoc)  = midval_martinb;
        estimPeeff(iLoc)    = midval_peeff;
        estimTeff2000(iLoc) = midval_teff2000;

        % Calculate z* and Teff 100 to 1000 m as a function of b 
        if (~isnan(midval_martinb))
            
%             % The following three methods produce similar results. I am
%             % going to use Method C as it requires less steps
%             
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
            
            [midval_teff1000,ci_teff1000,funvals_teff1000] = propagateErrorWithMC(funcTeff,funvals_martinb,'plot',false);
            [midval_zstar,ci_zstar,funvals_zstar] = propagateErrorWithMC(funcZstar,funvals_martinb,'plot',false);
        end

        if (nLocs < 10)
            fprintf('\nAnnual mean NPP is %4.3f and SST is %4.3f',qNppAnnualMean,qSstAnnualMean)
            fprintf('\nLocation %d',iLoc)
            fprintf('\nThe true b is %4.3f and the estimated b is %4.3f',trueMartinb(iLoc),estimMartinb(iLoc))
            fprintf('\nThe bounds estimated for b are %4.3f to %4.3f',ci_martinb(1),ci_martinb(2))
            fprintf('\nThe true PEeff is %4.3f and the estimated PEeff is %4.3f',truePeeff(iLoc),estimPeeff(iLoc))
            fprintf('\nThe bounds estimated for PEeff are %4.3f to %4.3f',ci_peeff(1),ci_peeff(2))
%             fprintf('\nThe estimated Teff 100 to 1000 m Method A is %4.3f, with bounds %4.3f to %4.3f',midval_teff1000_1,ci_teff1000_1(1),ci_teff1000_1(2))
%             fprintf('\nThe estimated Teff 100 to 1000 m Method B is %4.3f, with bounds %4.3f to %4.3f',midval_teff1000_2,ci_teff1000_2(1),ci_teff1000_2(2))
%             fprintf('\nThe estimated Teff 100 to 1000 m Method C is %4.3f, with bounds %4.3f to %4.3f\n',midval_teff1000_3,ci_teff1000_3(1),ci_teff1000_3(2))
            fprintf('\nThe estimated Teff 100 to 1000 m is %4.3f, with bounds %4.3f to %4.3f',midval_teff1000,ci_teff1000(1),ci_teff1000(2))
            fprintf('\nThe estimated z* is %4.0f, with bounds %4.0f to %4.0f\n',midval_zstar,ci_zstar(1),ci_zstar(2))
        end
        
        % Save arrays
        henson2012.martinb.ave(iLoc)       = midval_martinb;
        henson2012.martinb.stdevupp(iLoc)  = ci_martinb(2);
        henson2012.martinb.stdevlow(iLoc)  = ci_martinb(1);
        henson2012.martinb.max(iLoc)       = max(funvals_martinb);
        henson2012.martinb.min(iLoc)       = min(funvals_martinb);

        henson2012.teff1000.ave(iLoc)      = midval_teff1000;
        henson2012.teff1000.stdevupp(iLoc) = ci_teff1000(2);
        henson2012.teff1000.stdevlow(iLoc) = ci_teff1000(1);
        henson2012.teff1000.max(iLoc)      = max(funvals_teff1000);
        henson2012.teff1000.min(iLoc)      = min(funvals_teff1000);

        henson2012.zstar.ave(iLoc)         = midval_zstar;
        henson2012.zstar.stdevupp(iLoc)    = ci_zstar(2);
        henson2012.zstar.stdevlow(iLoc)    = ci_zstar(1);
        henson2012.zstar.max(iLoc)         = max(funvals_zstar);
        henson2012.zstar.min(iLoc)         = min(funvals_zstar);

    end % iLoc
    
    save(fullfile('.','data','processed','bcpmetrics_henson2012_local.mat'),'henson2012*')

end % isGlobalCalculation

% %% Contour plot of the variables controlling b
% 
% xtemp         = zeros(nr_ssta*nc_ssta,1); % for SST
% ynpp          = zeros(size(xtemp)); % for NPP
% zb            = zeros(size(xtemp)); % for b
% varflux2000   = zeros(size(xtemp));
% varexportprod = zeros(size(xtemp));
% 
% idx = 1;
% for iCol = 1:nc_ssta
%     for iRow = 1:nr_ssta
%         xtemp(idx) = qSstAnnualMean(iRow,iCol);
%         ynpp(idx) = qNppCarrAnnualMean(iRow,iCol);
%         zb(idx) = henson2012martinb(iRow,iCol);
%         varflux2000(idx) = henson2012ep2000(iRow,iCol);
%         varexportprod(idx) = henson2012ep100(iRow,iCol);
%         idx = idx + 1;
%     end
% end
% 
% [xtemp_sort, sortIdx] = sort(xtemp);
% ynpp_sort = ynpp(sortIdx);
% zb_sort = zb(sortIdx);
% 
% idxFirstNanTemp = find(isnan(xtemp_sort), 1, 'first');
% crop_xtemp_sort = xtemp_sort(1:idxFirstNanTemp-1);
% crop_ynpp_sort = ynpp_sort(1:idxFirstNanTemp-1);
% crop_ynpp_sort(isnan(crop_ynpp_sort)) = 0;
% crop_zb_sort = zb_sort(1:idxFirstNanTemp-1);
% crop_zb_sort(isnan(crop_zb_sort)) = 0;
% 
% mintemp = min(crop_xtemp_sort(:));
% maxtemp = max(crop_xtemp_sort(:));
% minnpp = min(crop_ynpp_sort(:));
% maxnpp = max(crop_ynpp_sort(:));
% 
% xv = linspace(mintemp, maxtemp, 150);
% yv = linspace(minnpp, maxnpp, 150);
% [Xm,Ym] = ndgrid(xv, yv);
% Zm = griddata(crop_xtemp_sort, crop_ynpp_sort, crop_zb_sort, Xm, Ym);
% 
% figure()
% set(gcf,'Units','Normalized','Position',[0.01 0.05 0.25 0.35],'Color','w') 
% contourf(Xm,Ym,Zm) % ,'ShowText','on'
% ylim([0 1000])
% xlim([-2.5 29])
% xh = xlabel('SST (ºC)');
% yh = ylabel('NPP (mg C m^{-2} d^{-1})');
% th = title('Henson et al. (2012) Martin b','FontSize',12);
% th.Position(2) = th.Position(2) + 30;
% xh.Position(2) = xh.Position(2) - 20;            
% yh.Position(1) = yh.Position(1) - 1;
% colorbar
% grid
% set(haxis(iSubplot),'FontSize',12)
% set(gcf,'PaperPositionMode','auto')
% print(gcf,fullfile(fullpathPlotsDir,'henson2012countour'),'-dpdf','-r0')  
% iFig = iFig + 1;
% 
% % [zb_sort2, sortIdx2] = sort(crop_zb_sort);
% % xtemp_sort2          = crop_xtemp_sort(sortIdx2);
% % zb_sort2             = crop_zb_sort(sortIdx2);
% % 
% % idxLastZeroZ       = find(zb_sort2 == 0, 1, 'last');
% % crop_xtemp_sort2   = xtemp_sort2(idxLastZeroZ+1:end);
% % crop_zb_sort2      = zb_sort2(idxLastZeroZ+1:end);
% % 
% % mdl = fitlm(crop_xtemp_sort2,crop_zb_sort2);
% % 
% % mymat = exportProduction(:,:);
% % integratedPOC = integrateAnnualPocAccordingToAreaOfTheGridCell(WOA13_lat,WOA13_lon,mymat);
% 
% end
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

% *************************************************************************

function [peeff] = Henson2012peeff(x)
    
    param.sst = x(1);
    param.sst_coeff1 = x(2);
    param.sst_coeff2 = x(3);

    % Particle export efficiency (Eq. S1 in the Supplementary Material of Henson et al. 2012)
    if (~isnan(param.sst))
        peeff = param.sst_coeff1*exp(-param.sst_coeff2*param.sst); % 0-1
    else
        peeff = 0;
    end
    
end % Henson2012peeff

% *************************************************************************

function [teff2000] = Henson2012teff2000(x)

    param.sst = x(1);
    param.sst_coeff1 = x(2);
    param.sst_coeff2 = x(3);
    param.npp = x(4);
    param.svi = x(5);

    % Export production
    if (~isnan(param.npp))
        ep = param.npp*Henson2012peeff(x); % mg m-2 d-1
    else
        ep = 0;
    end

    % The calculation of POC flux at 2000 m needs prd, rld and prr
    if (~isnan(param.svi))
        prd = 1e-3*(31*param.svi^2 + 49*param.svi + 7.8);
        rld = 1400*exp(-0.54*param.svi);
        prr = 1e-3*(2.6*param.svi^2 - 4.2*param.svi + 4.8);
    else
        prd = 0;
        rld = 0;
        prr = 0;
    end
    if (~isnan(param.npp))
        pocFlux2000 = param.npp*((prd*exp(-(2000-100)/rld))+prr); % mg m-2 d-1
    else
        pocFlux2000 = 0;
    end

    % Teff surface --> 2000 m
    if (ep > 0)
        teff2000 = pocFlux2000/ep;
    else
        teff2000 = 0;
    end

end % Henson2012teff2000

% *************************************************************************

function [martinb] = Henson2012martinb(x)

    teff2000 = Henson2012teff2000(x);

    % Martin's b, based on the power-law function F_2000 = EP x (2000/100)^b
    if (teff2000 > 0 && ~isnan(teff2000))
        martinb = -1.*(log(teff2000)/log(2000/100));
    else
        martinb = 0;
    end
    
end % Henson2012martinb

% *************************************************************************

end