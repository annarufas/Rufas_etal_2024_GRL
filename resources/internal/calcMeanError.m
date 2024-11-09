function [me,log_me] = calcMeanError(obs,pred)

% Also called mean difference (MD), mean error (ME), mean bias (MB) or just
% bias
% Greek symbol: delta

N = length(pred);

me = sum(obs-pred)/N; 

log_me = sum(log(obs)-log(pred))/N;

end