clc;
close all;
clear all

M = 500;
N = 1000;
A = randn(M,N);
% p = 20;
eps = 1e-3;

count = 1;
prange = 1:100;
errp = zeros(1,numel(prange));
distp = zeros(1,numel(prange));
for p = prange
    err = 0;
    dist = 0;
    for i = 1:100
        % Generate random support and input
        x0 = (2*(rand(N,1) > 0.5) - 1) .* rand(N,1);
        xs = randperm(N);
        x0(xs(p+1:N)) = 0;
        
        % Determine measurements
        b = A * x0;

        [x,xsupp] = omp(A,b,eps);
        
        % l2-norm error
        err = err + 0.01 * (norm(x-x0,2)/norm(x))^2;
        
        % Distance between original and output support
        num_xs = p;
        num_xsupp = numel(xsupp);
        xsupp_int = intersect(xs(1:p),xsupp);
        num_xsupp_int = numel(xsupp_int);
        dist_supp = 1 - (num_xsupp_int / max(num_xs,num_xsupp));
        dist = dist + 0.01 * dist_supp;
    end
    errp(count) = err;
    distp(count) = dist;
    count = count + 1;
%     fprintf('%d,%.4f\n',p,err);
end
figure, plot(prange,errp); 
figure, plot(prange,distp);