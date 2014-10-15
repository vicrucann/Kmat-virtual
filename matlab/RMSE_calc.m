function rmse = RMSE_calc(err)
n = size(err,2);
rmse = sqrt(sum(err.*err)/n);