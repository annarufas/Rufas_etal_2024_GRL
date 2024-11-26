function [fmartin,fexp] = solveMartinbAndZstar(isLogTransformed,fluxValues,fluxDepths)

% SOLVEMARTINBANDZSTAR This function solves for two key parameters related 
% to the attenuation of particulate organic carbon (POC) flux with depth: 
% Martin's b coefficient, which describes the rate of POC flux attenuation, 
% and the remineralisation length scale, z*, which represents the depth 
% scale at which 63% of POC flux ahs remineralised. The function uses a
% series of well-defined steps to ensure quality control of the fit.
% 
%   INPUT: 
%       isLogTransformed  - Boolean indicating whether to log-transform 
%                           depth and POC flux data before fitting b and z*
%       fluxValues        - POC flux values (mg C m-2 d-1)
%       fluxDepths        - POC flux cast depths (m)
%
%   OUTPUT:
%       fmartin - Martin's b coefficient, goodness of fit and standard
%                 deviation from the fit
%       fexp    - remineralisation length scale coefficient, goodness of 
%                 fit and standard deviation from the fit
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 23 Oct 2024 
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Define functions to solve for b and z*

if isLogTransformed

    % Linearise the power-law and exponential curves by applying log10
    
    fitMartinCurve = fittype('-b.*xb + a','problem',{'a'},'independent',{'xb'});
    fitExpCurve = fittype('-(1/zstar).*xz + a','problem',{'a'},'independent',{'xz'});

else 
    
    % Choose not to linearise the power-law and exponential curves to treat 
    % every data point as equally important.

    funcMartinCurve = @(b,f0,z0,z)     f0.*((z/z0).^(-b));      % Fz = F100 x (z/100)^b
    funcExpCurve    = @(zstar,f0,z0,z) f0.*exp(-(z-z0)./zstar); % Fz = F100 x e^(-(z-100)/z*)
    opts            = fitoptions('Method','NonlinearLeastSquares'); 

    fitMartinCurve = fittype(funcMartinCurve,'problem',{'f0','z0'},...
        'independent','z','coefficients','b','options',opts);
    fitExpCurve = fittype(funcExpCurve,'problem',{'f0','z0'},...
        'independent','z','coefficients','zstar','options',opts);
    
end

%% Initialise output arrays

fmartin = NaN(3,1); % [value, GOF, std of the fit]
fexp = NaN(3,1);

%% Start calculations

% The fit function cannot take NaN values, remove them
idxValid = ~isnan(fluxValues) & ~isnan(fluxDepths);
f = fluxValues(idxValid);
z = fluxDepths(idxValid);

% Proceed if there are data in the 1st position of the POC flux
% vector, there are at least 3 data points in the vertical
if (f(1) > 0 && size(f,1) >= 3) 

    f0 = f(1);
    z0 = z(1);

    if isLogTransformed

        a = log10(f0);
        xb = log10(z./z0);
        xz = z - z0;
        y = log10(f);

        [fm,gofm] = fit(xb,y,fitMartinCurve,'problem',{a},...
            'Lower',0.1,'Upper',4,'StartPoint',0.09);  
        [fe,gofe] = fit(xz,y,fitExpCurve,'problem',{a},...
            'Lower',50,'Upper',2000,'StartPoint',40); 

    else

        [fm,gofm] = fit(z,f,fitMartinCurve,'problem',{f0,z0},...
            'Lower',0.1,'Upper',4,'StartPoint',0.09);  
        [fe,gofe] = fit(z,f,fitExpCurve,'problem',{f0,z0},...
            'Lower',50,'Upper',2000,'StartPoint',40);

    end

    % Coefficient values
    fmartin(1) = fm.b;
    fmartin(2) = gofm.adjrsquare;
    fexp(1)    = fe.zstar;
    fexp(2)    = gofe.adjrsquare;

    % Calculate std from the 95% CI
    fmartin(3) = calculateStdFromCI(confint(fm), numel(fluxValues));
    fexp(3) = calculateStdFromCI(confint(fe), numel(fluxValues));
    
end

% Validate outputs
fmartin = validateFitOutput(fmartin, 0.10, 0.10, 4, 10);
fexp = validateFitOutput(fexp, 0.10, 50, 2000, 10);

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS USED IN THIS SCRIPT
% -------------------------------------------------------------------------

% *************************************************************************

function stdDev = calculateStdFromCI(coeffBounds,n)
    
    stdDev = (max(coeffBounds) - min(coeffBounds)) * sqrt(n) / 3.92;
    
end

% *************************************************************************

function validated = validateFitOutput(metrics,minGOF,minCoeff,maxCoeff,maxStdFactor)

    if isnan(metrics(2)) || isinf(metrics(2)) || metrics(2) < minGOF || ...
       metrics(1) < minCoeff || metrics(1) > maxCoeff || metrics(3) > maxStdFactor * metrics(1)
        validated = NaN(3,1); % return NaN array if validation fails
    else
        validated = metrics;
    end
    
end

% *************************************************************************

end % solveMartinbAndZstar