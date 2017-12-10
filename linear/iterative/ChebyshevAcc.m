function [Y, flag, iter, res] = ChebyshevAcc(R, C, rho, tol, maxIter, X0)
% ChebyshevAcc   Chebyshev acceleration.
%   X = ChebyshevAcc(R, C, RHO) attempts to find the fixed point of linear 
%   iteration X=R*X+C using Chebyshev acceleration. R is a N-by-N matrix 
%   whose eigenvalues are real numbers in (-1, 1). RHO is a scalar in (0,1)
%   which is larger than the spectrum radius of R.
%
%   X = ChebyshevAcc(R, C, RHO, TOL) specifies the tolerance of the method. 
%   If TOL is [] then ChebyshevAcc uses the default, 1e-6.
%
%   X = ChebyshevAcc(R, C, RHO, TOL, MAXITER) specifies the maximum number 
%   of iterations. If MAXITER is [] then ChebyshevAcc uses the default, 
%   max(N, 20).
%
%   X = ChebyshevAcc(R, C, RHO, TOL, MAXITER, X0) specifies the initial 
%   guess. If X0 is [] then ChebyshevAcc uses the default, an all zero 
%   vector.
%
%   [X, FLAG] = ChebyshevAcc(R, C, RHO, ...) also returns a convergence 
%   FLAG:
%    0 ChebyshevAcc converged to the desird tolerance TOL within MAXITER
%      iterations.
%    1 ChebyshevAcc iterated MAXITER times but did not converge.
%
%   [X, FLAG, ITER] = ChebyshevAcc(R, C, RHO, ...) also returns the 
%   iteration number at which X was computed: 0 <= ITER <= MAXITER.
%
%   [X, FLAG, ITER, RES] = ChebyshevAcc(R, C, RHO, ...) also returns a 
%   vector of the maximum of residual at each iteration. 
%
%   See also

%   Copyright 2017 Junshen Xu


CheckSquareMatrix(R);

if(~exist('maxIter', 'var') || isempty(maxIter))
    maxIter = 10;
end

if(~isreal(rho) || rho < 0 || rho > 1)
    error('rho should be a real number in (0, 1).');
end

if(~exist('X0','var') || isempty(X0))
    X0 = zeros(size(C));
end

CheckMultiplicationSize(R,X0, []);

if(~exist('tol', 'var') || isempty(tol))
    tol = 1e-6;
end

flag = 1;
res = zeros(maxIter, 1);

Y0 = X0;
Y1 = R * X0 + C;
res(1) = max(abs(Y1(:) - Y0(:)));
if(res(1) < tol)
    flag = 0;
    iter = 1;
    res = res(1);
    return
end

mu0 = 1;
mu1 = rho;
for iter = 2 : maxIter
    mu = 1 / (2 / (rho * mu1) - 1/ mu0);
    a = (2*mu / (rho * mu1));
    Y = R * Y1 * a - (mu / mu0) * Y0 + a * C;
    res(iter) = max(abs(Y(:) - Y1(:)));
    if(res(iter) < tol)
        flag = 0;
        break;
    end
    Y0 = Y1;
    Y1 = Y;
    mu0 = mu1;
    mu1 = mu;
end

res = res(1:iter);
