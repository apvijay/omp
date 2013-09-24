% Orthogonal matching pursuit (OMP)
% 
% Approximate the solution of P0: % min \|x\|_0 subject to Ax = b
% 
% Input  : Matrix A, vector b, error threshold eps
% Output : Vector x

function [x,supp] = omp(A,b,eps)

[M,N] = size(A);
% Find norm of columns of A
normcolA = zeros(1,N);
for j = 1:N
    normcolA(1,j) = norm(A(:,j),2);
end

% Initialization
% output
x = zeros(N,1);
% residual
r = b; % b - A*x;
% support
supp = [];

% figure; hold on;
% subplot(211), stem(x0,'bo');
% pause(1);

count = 0;
while norm(r,2) > eps
    count = count + 1;
    % Update support
    t = (A' * r).^2 ./ normcolA(:).^2;
    [~, ind] = max(t);
    supp = [supp ind];

    % Update residual
    As = A(:,supp);
    x1 = (As' * As) \ (As' * b);
    x(supp) = x1;
    r = b - A * x;
%     fprintf('%d,%.4f,%s\n',count,norm(r,2),num2str(supp));
    
%     subplot(212),stem(x,'r');
%     pause(0.7);
end
% (norm(x-x0,2)/norm(x))^2
