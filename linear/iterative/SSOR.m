function [X, flag, iter, res] = SSOR(A, B, omega, tol, maxIter, X0)
% SSOR   Solve linear system using symmetric successive over-relaxation.
%   X = SSOR(A, B, OMEGA) attempts to solve the linear system A*X=B for X,
%   where A is a N-by-N matrix and OMEGA is the relaxation factor. This
%   method is similar to SOR but will converge faster if A is symmetric.
%
%   X = SSOR(A, B, OMEGA, TOL) specifies the tolerance of the method. If 
%   TOL is [] then SSOR uses the default, 1e-6.
%
%   X = SSOR(A, B, OMEGA, TOL, MAXITER) specifies the maximum number of
%   iterations. If MAXITER is [] then SSOR uses the default, max(N, 20).
%
%   X = SSOR(A, B, OMEGA, TOL, MAXITER, X0) specifies the initial guess. If
%   X0 is [] then SSOR uses the default, an all zero vector.
%
%   [X, FLAG] = SSOR(A, B, OMEGA, ...) also returns a convergence FLAG:
%    0 SSOR converged to the desird tolerance TOL within MAXITER
%      iterations.
%    1 SSOR iterated MAXITER times but did not converge.
%    2 Jacobi Iterative matrix was ill-conditioned.
%
%   [X, FLAG, ITER] = SSOR(A, B, OMEGA, ...) also returns the iteration
%   number at which X was computed: 0 <= ITER <= MAXITER.
%
%   [X, FLAG, ITER, RES] = SSOR(A, B, OMEGA, ...) also returns a vector of
%   the maximum of residual at each iteration. 
%
%   See also SOR, JacIter, GSIter.

%   Copyright 2017 Junshen Xu

CheckSquareMatrix(A);
n = size(A, 1);
if(~exist('maxIter', 'var') || isempty(maxIter))
    maxIter = max(20,n);
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

if(omega < 0 || omega > 2)
    warning('omega out of [0, 2] may lead to divergence.')
end

if(min(abs(diag(A))) < eps)
    flag = 2;
    iter = 0;
    res = [];
    return
end

D = diag(diag(A));
invDL = Inverse(D + omega * tril(A, -1));
invDU = Inverse(D + omega * triu(A, 1));
BSOR1 = invDL * ((1-omega)*D - omega * triu(A, 1)); 
BSOR2 = invDU * ((1-omega)*D - omega * tril(A, -1));
FSOR1 = omega * invDL * B;
FSOR2 = omega * invDU * B;

flag = 1;
res = zeros(maxIter, 1);
for iter = 1 : maxIter
    X = BSOR1 * X + FSOR1;
    X = BSOR2 * X + FSOR2;
    res(iter) = max(max(abs(A*X - B)));
    if(res(iter) < tol)
        flag = 0;
        break;
    end
end

res = res(1:iter);
