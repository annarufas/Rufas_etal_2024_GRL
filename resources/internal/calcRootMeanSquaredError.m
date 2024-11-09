function [rmse,log_rmse] = calcRootMeanSquaredError(obs,pred)

N = length(obs);

squaredErrors = (obs-pred).^2;
sumOfSquaredErrors = sum(squaredErrors);
rmse = sqrt(sumOfSquaredErrors/N);

squaredLogErrors = (log(obs)-log(pred)).^2; 
sumOfSquaredLogErrors = sum(squaredLogErrors);
log_rmse = sqrt(sumOfSquaredLogErrors/N);

end