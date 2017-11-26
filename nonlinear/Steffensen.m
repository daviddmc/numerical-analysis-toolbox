function x = Steffensen(phi, x0, maxIter, tolerance)
% Steffensen acceleration for fixed-point iteration
% input
% phi : iteration function
% x0 : initail value of x, default: 0
% maxIter : max number of iterations, default: 10
% tolerance : the tolerance of the method, default: 1e-6
% output
% x

if(~exist('x0','var') || isempty(x0))
    x0 = 0;
end

if(~exist('maxIter','var') || isempty(maxIter))
    maxIter = 10;
end

if(~exist('tolerance','var') || isempty(tolerance))
    tolerance = 1e-6;
end

x = x0;

for k = 1 : maxIter
    y = phi(x);
    z = phi(y);
    x_new = x - (y-x)^2/(z-2*y+x);
    if(abs(x_new - x)<tolerance)
        x = x_new;
        break;
    end
    x = x_new;
end


