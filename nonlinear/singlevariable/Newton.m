function [x, flag, iter, xs] = Newton(fun, grad, x0, tol, maxIter)
% Newton   Find zero of single-variable function using Newton method.
%   X = Newton(FUN, GRAD, X0) returns a zero of FUN using Newton method.
%   FUN should be a function handle taking a scalar as input. GRAD is 
%   the derivative of FUN. X0 is the initial point which is a scalar. 
%
%   X = Newton(FUN, GRAD, X0, TOL) specifies the tolerance of the 
%   method. If TOL is [] then Newton uses the default, 1e-6. 
%   
%   X = Newton(FUN, GRAD, X0, TOL, MAXITER) specifies the maximum 
%   number of iterations. If MAXITER is [] then Newton uses the 
%   default, 10.
%
%   [X, FLAG] = Newton(FUN, GRAD, X0, ...) also returns a convergence 
%   FLAG:
%    0 Newton converged to the desird tolerance TOL within MAXITER
%      iterations.
%    1 Newton iterated MAXITER times but did not converge.
%
%   [X, FLAG, ITER] = Newton(FUN, GRAD, X0, ...) also returns the 
%   iteration number at which X was computed: 0 <= ITER <= MAXITER.
%
%   [X, FLAG, ITER, XS] = Newton(FUN, GRAD, X0 ...) also returns a 
%   vector of iteration result at each step, length(XS) = ITER, 
%   XS(end) = X;
%
%   See also

%   Copyright 2017 Junshen Xu

if(~exist('x0','var') || isempty(x0))
    x0 = 0;
end

if(~exist('maxIter','var') || isempty(maxIter))
    maxIter = 10;
end

if(~exist('tol','var') || isempty(tol))
    tol = 1e-6;
end



%if(grad(x0) == 0)
%    x0 = x0 + 200 * rand * tol;
%end

flag = 1;
fx = fun(x0);
dfx = grad(x0);
xs = zeros(maxIter, 1);
x = x0;
for iter = 1 : maxIter
    if abs(fx) < tol
        flag = 0;
        xs(iter) = x;
        break;
    elseif abs(dfx) < tol
        x0 = x0 + 200 * rand * tol;
    end
    
    delta = fx / dfx;
    x = x0 - delta;
    x1 = x - delta;
    fx1 = fun(x1);
    fx = fun(x);
    while abs(fx1) < abs(fx)
        x = x1;
        fx = fx1;
        x1 = x - delta;
        fx1 = fun(x1);
    end
    
    dfx = grad(x);
    xs(iter) = x;
    x0 = x;
end

xs = xs(1:iter);

%{
if ~flag
    for k = 1 : maxIter2
        fx = f(x0);
        dfx = df(x0);
        delta = fx * log(abs(fx)) / dfx / (log(abs(fx)) - log(abs(dfx)));
        if(isnan(delta))
            flag = 1;
        else
            x = x0 - delta;
        end
        if(abs(delta) < tolerance || (abs(fx) < tolerance && abs(dfx) < tolerance))
            flag = 1;
        end
        if(flag)
            break;
        end
        x0 = x;
    end
end
%}