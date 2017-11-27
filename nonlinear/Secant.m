function x = Secant( f, x0, x1, maxIter, tolerance)
% Secant method

if(~exist('x0','var') || isempty(x0))
    x0 = 0;
end

if(~exist('x1','var') || isempty(x1))
    x1 = 1;
end

if(~exist('maxIter','var') || isempty(maxIter))
    maxIter = 10;
end

if(~exist('tolerance','var') || isempty(tolerance))
    tolerance = 1e-6;
end


for k = 1 : maxIter
    y1 = f(x1);
    y0 = f(x0);
    x = x1 - y1 * (x1 - x0) / (y1 - y0);
    if(abs(x - x1) < tolerance)
        break
    end
    x0 = x1;
    x1 = x;
end