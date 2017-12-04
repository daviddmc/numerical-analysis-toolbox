function [X, flag, iter, res] = GSIter(A, B, tol, maxIter, X0)
% GSIter   Solve linear system using Gauss-Seidel iteration.
%   X = GSIter(A, B) attempts to solve the linear system A*X=B for X,
%   where A is a N-by-N matrix.
%
%   X = GSIter(A, B, TOL) specifies the tolerance of the method. If TOL is
%   [] then GSIter uses the default, 1e-6.
%
%   X = GSIter(A, B, TOL, MAXITER) specifies the maximum number of
%   iterations. If MAXITER is [] then GSIter uses the default, min(N, 20).
%
%   X = GSIter(A, B, TOL, MAXITER, X0) specifies the initial guess. If X0
%   is [] then GSIter uses the default, an all zero vector.
%
%   [X, FLAG] = GSIter(A, B, ...) also returns a convergence FLAG:
%    0 GSIter converged to the desird tolerance TOL within MAXITER
%      iterations.
%    1 GSIter iterated MAXITER times but did not converge.
%    2 Gauss-Seidel Iterative matrix was ill-conditioned.
%
%   [X, FLAG, ITER] = GSIter(A, B, ...) also returns the iteration number
%   at which X was computed: 0 <= ITER <= MAXITER.
%
%   [X, FLAG, ITER, RES] = GSIter(A, B, ...) also returns a vector of the
%   maximum of residual at each iteration. 
%
%   See also

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

if(min(abs(diag(A))) < eps)
    flag = 2;
    X = X0;
    iter = 0;
    res = [];
    return
end
invDL = Inverse(tril(A, 0));

BG = eye(n) - invDL * A;
FG = invDL * B;

res = zeros(maxIter, 1);
flag = 1;
for iter = 1 : maxIter
    X = BG * X + FG;
    res(iter) = max(max(abs(A*X - B)));
    if(res(iter) < tol)
        flag = 0;
        break;
    end
end

res = res(1:iter);
