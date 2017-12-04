function [m, u, flag, iter] = RayleighQuotientIter(A, tol, maxIter, u0)
% RayleighQuotientIter     Rayleigh Quotient iteration method.
%   M = RayleighQuotientIter(A) tries to find the dominant eigenvalue of a 
%   (conjugate) symmetric matrix A.
%
%   M = RayleighQuotientIter(A, TOL) specifies the tolerance of the method. 
%   If TOL is [] then RayleighQuotientIter uses the default, 1e-6.
%
%   M = RayleighQuotientIter(A, TOL, MAXITER) specifies the maximum number 
%   of iterations. If MAXITER is [] then RayleighQuotientIter uses the 
%   default, max(N, 20), where N is the size of A.
%
%   M = RayleighQuotientIter(A, TOL, MAXITER, U0) specifies the initial 
%   guess of the eigenvector of the dominant eigenvalue. If U0 is [] then 
%   RayleighQuotientIter uses the default, an all one vector.
%
%   [M, U] = RayleighQuotientIter(A, ...) also returns the corresponding 
%   eigenvector.
%
%   [M, U, FLAG] = RayleighQuotientIter(A, ...) also returns a convergence 
%   FLAG:
%    0 RayleighQuotientIter converged to the desird tolerance TOL within 
%      MAXITER iterations.
%    1 RayleighQuotientIter iterated MAXITER times but did not converge.
%
%   [M, U, FLAG, ITER] = RayleighQuotientIter(A, ...) also returns the 
%   iteration number at which M was computed: 0 <= ITER <= MAXITER.
%
%   See also InverseIter.

%   Copyright 2017 Junshen Xu

CheckSquareMatrix(A);
n = size(A);
if(~exist('maxIter', 'var') || isempty(maxIter))
    maxIter = max(20, n);
end

if(~exist('u0','var') || isempty(u0))
    u0 = ones(size(A, 1), 1);
end

u = u0 / Norm(u0);

CheckMultiplicationSize(A, u0, []);

if(~exist('tol', 'var') || isempty(tol))
    tol = 1e-6;
end

flag = 1;
for iter = 1 : maxIter
    m = real(u' * A * u);
    u = GaussianElimination(A - m * eye(n), u, 'complete');
    u = u / Norm(u);
    if(Norm(A*u - m * u) < tol)
        flag = 0;
        break
    end
end


end

