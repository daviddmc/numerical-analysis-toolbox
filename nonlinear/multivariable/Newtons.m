function [x, flag, iter] = Newtons(fun, jac, x0, tol, maxIter)
% Newtons   Solving nonlinear equations with Newton method.
%   X = Newtons(FUN, JAC, X0) starts at the column vector X0 and tries to 
%   solve the equations in FUN.  FUN accepts input X and returns a column 
%   vector of equation values FUN evaluated at X. JAC accepts input X and 
%   returns the Jacobian matrix at X.
%
%   X = Newtons(FUN, JAC, X0, TOL) specifies the tolerance. If TOL is [],
%   Newtons uses the default, 1e-6.
%
%   X = Newtons(FUN, JAC, X0, TOL, MAXITER) specifies the max number of 
%   iterations. If MAXITER is [], Newtons uses the default, 10. 
%   
%   [X, FLAG] = Newtons(FUN, JAC, X0, ...) also returns a convergence FLAG:
%
%     0 Newtons converged to the desird tolerance TOL within MAXITER
%       iterations.
%     1 Newtons iterated MAXITER times but did not converge.
%
%   [X, FLAG, ITER] = Newtons(FUN,JAC, X0, ...) also returns the iteration 
%   number at which X was computed: 0 <= ITER <= MAXITER.
%
%   See also Broyden, Continuation.

%   Copyright 2017 Junshen Xu

if ~exist('tol','var') || isempty(tol)
    tol = 1e-6;
end

if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter = 10;
end

x = x0;
flag = 1;
for iter = 1:maxIter
    f = fun(x);
    J = jac(x);
    y = GaussianElimination(J, -f);
    x = x + y;
    if Norm(y) < tol
        flag = 0;
        break
    end
end

