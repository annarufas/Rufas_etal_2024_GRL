function mrpe = calcMeanRelativePercentageError(obs,pred)

N = length(pred);
relativePercentageError = 100.*((obs-pred)./obs);
sumOfRelativePercentageError = sum(relativePercentageError);
mrpe = sumOfRelativePercentageError/N;

end
