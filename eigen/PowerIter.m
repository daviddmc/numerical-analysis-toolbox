function [m, u, flag, iter] = PowerIter( A, tol, maxIter, u0)
% PowerIter     Power iteration method.
%   M = PowerIter(A) tries to find the dominant eigenvalue of matrix A.
%
%   M = PowerIter(A, TOL) specifies the tolerance of the method. If 
%   TOL is [] then PowerIter uses the default, 1e-6.
%
%   M = PowerIter(A, TOL, MAXITER) specifies the maximum number of
%   iterations. If MAXITER is [] then PowerIter uses the default, 
%   max(N, 20).
%
%   M = PowerIter(A, TOL, MAXITER, U0) specifies the initial guess of the 
%   eigenvector of the dominant eigenvalue. If U0 is [] then PowerIter uses 
%   the default, an all one vector.
%
%   [M, U] = PowerIter(A, ...) also returns the corresponding eigenvector.
%
%   [M, U, FLAG] = PowerIter(A, ...) also returns a convergence FLAG:
%    0 PowerIter converged to the desird tolerance TOL within MAXITER
%      iterations.
%    1 PowerIter iterated MAXITER times but did not converge.
%
%   [M, U, FLAG, ITER] = PowerIter(A, ...) also returns the iteration
%   number at which M was computed: 0 <= ITER <= MAXITER.
%
%   See also InverseIter.

%   Copyright 2017 Junshen Xu

CheckSquareMatrix(A);
n = size(A, 1);
if(~exist('maxIter', 'var') || isempty(maxIter))
    maxIter = max(20, n);
end

if(~exist('u0','var') || isempty(u0))
    u0 = ones(size(A, 1), 1);
end

CheckMultiplicationSize(A, u0, []);

if(~exist('tol', 'var') || isempty(tol))
    tol = 1e-6;
end

u = u0;
flag = 1;
for iter = 1 : maxIter
    v = A * u;
    [~, idx] = max(abs(v));
    m = v(idx);
    if(max(abs(v - m*u)) < tol)
        flag = 0;
        break
    end
    u = v / m;
end



