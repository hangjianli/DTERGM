function res = admm_gfl(Y, lam)

[n , ~] = size(Y);
option.weights = ones(n-1);
option.verbose = 1;
option.max_time = 1000;

res = gflasso(Y, lam, option);

end

