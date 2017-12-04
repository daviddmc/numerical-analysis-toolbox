function [X, flag, iter, res] = SOR(A, B, omega, tol, maxIter, X0)
% SOR   Solve linear system using successive over-relaxation.
%   X = SOR(A, B, OMEGA) attempts to solve the linear system A*X=B for X,
%   where A is a N-by-N matrix and OMEGA is the relaxation factor.
%
%   X = SOR(A, B, OMEGA, TOL) specifies the tolerance of the method. If TOL 
%   is [] then SOR uses the default, 1e-6.
%
%   X = SOR(A, B, OMEGA, TOL, MAXITER) specifies the maximum number of
%   iterations. If MAXITER is [] then SOR uses the default, min(N, 20).
%
%   X = SOR(A, B, OMEGA, TOL, MAXITER, X0) specifies the initial guess. If
%   X0 is [] then SOR uses the default, an all zero vector.
%
%   [X, FLAG] = SOR(A, B, OMEGA, ...) also returns a convergence FLAG:
%    0 SOR converged to the desird tolerance TOL within MAXITER
%      iterations.
%    1 SOR iterated MAXITER times but did not converge.
%    2 Jacobi Iterative matrix was ill-conditioned.
%
%   [X, FLAG, ITER] = SOR(A, B, OMEGA, ...) also returns the iteration
%   number at which X was computed: 0 <= ITER <= MAXITER.
%
%   [X, FLAG, ITER, RES] = SOR(A, B, OMEGA, ...) also returns a vector of
%   the maximum of residual at each iteration. 
%
%   See also

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

if(omega < 0 || omega > 2)
    warning('omega out of [0, 2] may lead to divergence.')
end

CheckMultiplicationSize(A,X0, B);

if(~exist('tol', 'var') || isempty(tol))
    tol = 1e-6;
end
tol = max(tol, tol * max(abs(B(:))));

if(min(abs(diag(A))) < eps)
    flag = 2;
    iter = 0;
    res = [];
    return
end

D = diag(diag(A));
invDL = Inverse(D + omega * tril(A, -1));
BSOR = invDL * ((1-omega)*D - omega * triu(A, 1)); 
FSOR = omega * invDL * B;

flag = 1;
res = zeros(maxIter,1);
for iter = 1 : maxIter
    X = BSOR * X + FSOR;
    res(iter) = max(max(abs(A*X - B)));
    if(res(iter) < tol)
        flag = 0;
        break;
    end
end

res = res(1 : iter);
