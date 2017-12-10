function [x, flag, iter, xs]  = Muller(fun, x0, x1, x2, tol, maxIter)
% Muller   Find zero of single-variable function using Muller method.
%   X = Muller(FUN, X0, X1, X2) returns a zero of FUN using Muller 
%   method. FUN should be a function handle taking a scalar as input.
%   X0, X1, X3 are the initial points which are scalar. 
%
%   X = Muller(FUN, X0, X1, X2 TOL) specifies the tolerance of the 
%   method. If TOL is [] then Muller uses the default, 1e-6. 
%   
%   X = Muller(FUN, X0, X1, X2, TOL, MAXITER) specifies the maximum 
%   number of iterations. If MAXITER is [] then Muller uses the 
%   default, 10.
%
%   [X, FLAG] = Muller(FUN, X0, X1, X2 ...) also returns a 
%   convergence FLAG:
%    0 Muller converged to the desird tolerance TOL within MAXITER
%      iterations.
%    1 Muller iterated MAXITER times but did not converge.
%
%   [X, FLAG, ITER] = Muller(FUN, X0, X1, X2 ...) also returns the 
%   iteration number at which X was computed: 0 <= ITER <= MAXITER.
%
%   [X, FLAG, ITER, XS] = Muller(FUN, X0, X1, X2 ...) also returns 
%   a vector of iteration result at each step, length(XS) = ITER, 
%   XS(end) = X;
%
%   See also

%   Copyright 2017 Junshen Xu


if(~exist('maxIter','var') || isempty(maxIter))
    maxIter = 10;
end

if(~exist('tol','var') || isempty(tol))
    tol = 1e-6;
end

%x2 = x0;
%x1 = x2 + 100 * tolerance;
%x0 = x2 - 100 * tolerance;
x = x2;

flag = 1;
flagRecord = nargout > 3;
if flagRecord
    xs = zeros(maxIter,1);
end

for iter = 1 : maxIter
    h1 = x1 - x0;
    h2 = x2 - x1;
    f0 = fun(x0);
    f1 = fun(x1);
    f2 = fun(x2);
    delta1 = (f1 - f0) / h1;
    delta2 = (f2 - f1) / h2;
    d = (delta2 - delta1) / (h2 + h1);
    
    b = delta2 + h2 * d;
    D = sqrt(b^2 - 4*f2*d);
    
    if(abs(b-D) < abs(b+D))
        E = b + D;
    else
        E = b - D;
    end
    
    h = -2 * f2 / E;
    if(isnan(h) || isinf(h))
        flag = 0;
    else
        x = x2 + h;
    end
    
    if(abs(h) < tol)
        flag = 0;
    end

    if flagRecord
        xs(iter) = x;
    end
    
    if(~flag)
        break;
    end
    
    x0 = x1;
    x1 = x2;
    x2 = x;
end

if flagRecord
    xs = xs(1:iter);
end