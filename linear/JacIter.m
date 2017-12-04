function [X, flag, iter, res] = JacIter(A, B, tol, maxIter, X0)
% JacIter   Solve linear system using Jacobi iteration.
%   X = JacIter(A, B) attempts to solve the linear system A*X=B for X,
%   where A is a N-by-N matrix.
%
%   X = JacIter(A, B, TOL) specifies the tolerance of the method. If TOL is
%   [] then JacIter uses the default, 1e-6.
%
%   X = JacIter(A, B, TOL, MAXITER) specifies the maximum number of
%   iterations. If MAXITER is [] then JacIter uses the default, min(N, 20).
%
%   X = JacIter(A, B, TOL, MAXITER, X0) specifies the initial guess. If X0
%   is [] then JacIter uses the default, an all zero vector.
%
%   [X, FLAG] = JacIter(A, B, ...) also returns a convergence FLAG:
%    0 JacIter converged to the desird tolerance TOL within MAXITER
%      iterations.
%    1 JacIter iterated MAXITER times but did not converge.
%    2 Jacobi Iterative matrix was ill-conditioned.
%
%   [X, FLAG, ITER] = JacIter(A, B, ...) also returns the iteration number
%   at which X was computed: 0 <= ITER <= MAXITER.
%
%   [X, FLAG, ITER, RES] = JacIter(A, B, ...) also returns a vector of the
%   maximum of residual at each iteration. 
%
%   See also GSIter, SOR, SSOR

%   Copyright 2017 Junshen Xu

CheckSquareMatrix(A);
n = size(A, 1);

if(~exist('maxIter', 'var') || isempty(maxIter))
    maxIter = max(20, n);
end

if(~exist('X0','var') || isempty(X0))
    X0 = zeros(size(B));
end
X = X0;

CheckMultiplicationSize(A,X0, B);

if(~exist('tol', 'var') || isempty(tol))
    tol = 1e-6;
end

tol = max(tol, tol * max(abs(B(:))));

D = diag(A);
if(min(abs(D)) < eps)
    flag = 2;
    iter = 0;
    res = [];
    return
end
flag = 1;
BJ = eye(n);
FJ = B;
for ii = 1 : n
    BJ(ii, :) = BJ(ii, :) - A(ii, :) / D(ii);
    FJ(ii, :) = B(ii, :) / D(ii);
end

res = zeros(maxIter, 1);
for iter = 1 : maxIter
    X = BJ * X + FJ;
    res(iter) = max(max(abs(A*X - B)));
    if(res(iter) < tol)
        flag = 0;
        break;
    end
end

res = res(1:iter);

