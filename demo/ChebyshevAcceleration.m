% This script shows the effect of Chebyshev acceleration in iterative 
% method for solving linear system. For more infomation, see section 6.5.6 
% in Applied Numerical Linear Algebra, James W. Demmel.
%
% Copyright 2017 Junshen Xu

%% setup
N = 10; % matrix size
A = 3 * eye(N) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
b = zeros(N, 1);
x0 = ones(N, 1); 
% jacobi iterative matrix
D = diag(diag(A));
BJ = eye(N) - D\A;
FJ = zeros(N, 1);
% calcualte rho
rho = max(abs(SymEigen(BJ)));
tol = 1e-6;
maxIter = ceil(-log(tol) / rho * 2);

%% iteration
[~, flagJac, iterJac, resJac] = JacIter(A, b, tol, maxIter, x0);
[~, flagChe, iterChe, resChe] = ChebyshevAcc(BJ, FJ, rho, tol, maxIter, x0);

%% show results
semilogy(1:iterJac, resJac, 1:iterChe, resChe);
hold on
plot(1:max([iterJac,iterChe]), 1e-6*ones(max([iterJac,iterChe]),1))
xlabel('\# of iteration', 'Interpreter','latex')
ylabel('$error$', 'Interpreter','latex');
legend('Jcobi','Chebyshev','location','best');
title('Chebyshev Acceleration','Interpreter','latex');
grid on