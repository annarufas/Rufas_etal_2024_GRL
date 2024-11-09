function mape = calcMeanAbsolutePercentageError(obs,pred)

% Also called mean absolute percentage difference (MAPD) or mean absolute
% percentage error (MAPE)
N = length(pred);
absolutePercentageError = 100.*(abs(obs-pred)./obs);
sumOfAbsolutePercentageError = sum(absolutePercentageError);
mape = sumOfAbsolutePercentageError/N;

end
