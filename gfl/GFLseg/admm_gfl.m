function res = admm_gfl(Y, lam, tol, maxit)

[n , ~] = size(Y);
option.weights = ones(n-1);
option.verbose = 1;
option.max_time = 1000;
option.tol = tol;
option.maxit = maxit;

res = gflasso(Y, lam, option);

end

