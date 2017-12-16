function [m, u, flag, iter] = PowerIter( A, tol, maxIter, u0, symString)
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
%   M = PowerIter(A, TOL, MAXITER, U0, 'symmetric') assumes that A is
%   symmetric(Hermitian) matrix.
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

flag = 1;
u = u0;
m0 = 0;
m1 = 0;
if(exist('symString','var') && strncmpi(symString,'sym',1))
    u = u / Norm(u);
    for iter = 1:maxIter
        v = A*u;
        m = u'*v;
        vNorm = Norm(v);
        if(Norm(v/vNorm - u) < tol)
            flag = 0;
            if iter>=4
                m = m0 - (m1 - m0)^2 / (m -2*m1 + m0);
            end
            break;
        end
        u = v / vNorm;
        m0 = m1;
        m1 = m;
    end
else
    [~, idx] = max(abs(u));
    u = u / u(idx);
    for iter = 1 : maxIter
        v = A * u;
        [~, idx] = max(abs(v));
        m = v(idx);
        if(max(abs(v/m - u)) < tol)
            flag = 0;
            if iter >=4
                m = m0 - (m1 - m0)^2 / (m -2*m1 + m0);
            end
            break
        end
        u = v / m;
        m0 = m1;
        m1 = m;
    end
end


