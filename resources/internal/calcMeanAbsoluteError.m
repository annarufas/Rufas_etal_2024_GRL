function [mae,log_mae] = calcMeanAbsoluteError(obs,pred)

% Also called mean absolute deviation or difference (MAD) or error (MAE)

N = length(pred);

mae = sum(abs(obs-pred))/N; 

log_mae = sum(abs(log(obs)-log(pred)))/N;

end