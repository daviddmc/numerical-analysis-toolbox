function [x, flag] = Muller( f, x0, maxIter, tolerance)
% Newton method

flag = 0;

if(~exist('x0','var') || isempty(x0))
    x0 = 0;
end

if(~exist('maxIter','var') || isempty(maxIter))
    maxIter = 10;
end

if(~exist('tolerance','var') || isempty(tolerance))
    tolerance = 1e-6;
end

x2 = x0;
x1 = x2 + 100 * tolerance;
x0 = x2 - 100 * tolerance;
x = x2;

for k = 1 : maxIter
    h1 = x1 - x0;
    h2 = x2 - x1;
    f0 = f(x0);
    f1 = f(x1);
    f2 = f(x2);
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
        flag = 1;
    else
        x = x2 + h;
    end
    
    if(abs(h) < tolerance)
        flag = 1;
    end
    
    if(flag)
        break;
    end
    
    x0 = x1;
    x1 = x2;
    x2 = x;
end