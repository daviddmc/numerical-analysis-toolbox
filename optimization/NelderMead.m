function [x, y] = NelderMead(f, x0, rad, maxIter, epsilon, delta)
%NELDERMEAD Summary of this function goes here
%   Detailed explanation goes here
n = length(x0);
fac = factorial(n);
x(:, 1) = x0;
x(:, 2:n+1) = x0*ones(1,n) + rad*eye(n);
for j = n+1 : -1 : 1
    y(j) = f(x(:, j));
end
[y, r] = sort(y);
x = x(:, r);
for i = 1 : maxIter
    xbar = mean(x(:, 1:n), 2);
    xh = x(:, n+1);
    xr = 2 * xbar - xh;
    yr = f(xr);
    if yr < y(n)
        if yr < y(1)
            xe = 3*xbar - 2*xh;
            ye = f(xe);
            if ye < yr
                x(:, n+1) = xe;
                y(n+1) = ye;
            else
                x(:, n+1) = xr;
                y(n+1) = yr;
            end
        else
            x(:, n+1) = xr;
            y(n+1) = yr;
        end
    else
        if yr < y(n+1)
            xoc = 1.5 * xbar - 0.5 * xh;
            yoc = f(xoc);
            if yoc < yr
                x(:,n+1) = xoc;
                y(n+1) = yoc;
            else
                for j = 2 : n+1
                    x(:,j) = 0.5*x(:,1) + 0.5*x(:,j);
                    y(j) = f(x(:,j));
                end
            end
        else
            xic = 0.5*xbar + 0.5*xh;
            yic = f(xic);
            if yic < y(n+1)
                x(:,n+1) = xic;
                y(n+1) = yic;
            else
                for j = 2:n+1
                    x(:,j) = 0.5 *x(:,1) + 0.5*x(:,j);
                    y(j) = f(x(:,j));
                end
            end
        end
    end
    [y, r] = sort(y);
    x = x(:, r);
    if(y(end) - y(1) < epsilon)
        xo = bsxfun(@minus, x(:,2:end), x(:,1));
        if(det(xo) < delta * fac)
            x = x(:, 1);
            y = y(1);
            break;
        end
    end
end

