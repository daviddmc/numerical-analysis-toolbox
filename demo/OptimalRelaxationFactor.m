% This example shows the effect of relaxation factor on SOR iteration.
% Consider the Poisson problem in a N-by-N gird, which can be formulated as
% the following N^2-by-N^2 linear system P*x=b. P is the following matrix
%
%   [  A  -I        ]
%   [ -I   A  \\    ]
%   [     \\  \\  -I]
%   [         -I   A]
%
% where A is N-by-N tridiagonal matrix with elements in main diagonal to
% be 4 and elements in super- and sub- diagonal to be -1. I is N-by-N
% identity matrix.
%
% According to Young's estimation, the optimal relaxation factor of P is  
% 2 / (1 + sin(pi / (n+1))).
%
% This script evaluates the speed of convergence of 3 different methods:
%   1. Jacobi iteration
%   2. Gauss-Seidel iteration (SOR with omega = 1)
%   3. SOR iteration (with omega = omega_optimal)

%   Copyright 2017 Junshen Xu

%% parameters
n = 9; % size of grid
epsilon = 1e-6; % tolerance of error

%% setup
maxIter = ceil(-log10(epsilon) * log(10) * (n+1)^2 / pi^2); % max # of iteration
% calculate maxtrix A
block = 4 * eye(n) + diag(-ones(n-1,1),1) + diag(-ones(n-1,1),-1);
blocks = cell(n, 1);
blocks(:) = {block};
A = blkdiag(blocks{:}) - diag(ones(n^2-n,1),n) -diag(ones(n^2-n,1),-n);
% set omega to the optimal value
omega = 2 / (1 + sin(pi / (n+1)));
% init
x0 = ones(n^2, 1);
b = zeros(n^2, 1);

%% iteration
[~, ~, iterJac, resJac] = JacIter(A, b, epsilon, maxIter, x0);
[~, ~, iterGS, resGS]  = GSIter(A, b, epsilon, maxIter, x0);
[~, ~, iterSOR, resSOR] = SOR(A, b, omega, epsilon, maxIter, x0);

%% show results
semilogy(1:iterJac, resJac, 1:iterSOR, resSOR, 1:iterGS, resGS);
xlabel('\# of iteration', 'Interpreter','latex')
ylabel('$||x-x^*||_\infty$', 'Interpreter','latex');
legend('Jcobi','SOR','Gauss-Seidel','location','best');
text(iterJac,resJac(end),['(' num2str(iterJac) ', ' num2str(resJac(end)) ')'],...
    'HorizontalAlignment', 'right');
text(iterGS,resGS(end),['(' num2str(iterGS) ', ' num2str(resGS(end)) ')'],...
     'HorizontalAlignment', 'center');
text(iterSOR,resSOR(end),['(' num2str(iterSOR) ', ' num2str(resSOR(end)) ')']);


