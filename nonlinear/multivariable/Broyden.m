function [x, flag, iter] = Broyden(fun, x0, tol, maxIter)
% Broyden   Broyden's rank 1 quasi-Newton method.
%   X = Broyden(fun, X0) starts at the column vector X0 and tries to solve 
%   the equations in FUN.  FUN accepts input X and returns a column vector 
%   of equation values FUN evaluated at X.
%
%   X = Broyden(fun, X0, TOL) specifies the tolerance. If TOL is [],
%   Broyden uses the default, 1e-6.
%
%   X = Broyden(fun, X0, TOL, MAXITER) specifies the max number of 
%   iterations. If MAXITER is [], Broyden uses the default, 10. 
%   
%   [X, FLAG] = Broyden(FUN, X0, ...) also returns a convergence FLAG:
%
%     0 Broyden converged to the desird tolerance TOL within MAXITER
%       iterations.
%     1 Broyden iterated MAXITER times but did not converge.
%
%   [X, FLAG, ITER] = Broyden(FUN, X0, ...) also returns the iteration 
%   number at which X was computed: 0 <= ITER <= MAXITER.
%
%   See also Newtons, Continuation.

%   Copyright 2017 Junshen Xu

if ~exist('tol','var') || isempty(tol)
    tol = 1e-6;
end

if ~exist('matIter', 'var') || isempty(maxIter)
    maxIter = 10;
end

n = length(x0);
fx = fun(x0);
Binv = eye(n);
x = x0;
flag = 1;
for iter = 1 : maxIter
    if Norm(fx) < tol
        flag = 0;
        break
    end
    s = -Binv*fx;
    xnew = x + s;
    fnew = fun(xnew);
    y = fnew - fx;
    Binvy = Binv*y;
    Binv = Binv + (s - Binvy)/(s'*Binvy)*(s'*Binv);
    x = xnew;
    fx = fnew;
end

end

