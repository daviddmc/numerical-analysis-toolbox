function a = Backtracking(fun, grad, x, p, a0, rho, c)
%BACKTRACKING Summary of this function goes here
%   Detailed explanation goes here

% Armijo
% rho in (0,1)
% c in (0,1)

fx = fun(x);
a = a0;
fxp = f(x + a*p);
b = grad'*p*c;
while fxp > fx + b*a
    a = rho * a;
    fxp = f(x + a*p);
end

function Wolfe(fun, grad, x, p, a0, rho, c1, c2)
% 0 < c1 < c2 < 1
strong = 1;
a = a0;
fx = fun(x);
gx = grad(x);
b = gx'*p;
c1 = c1 * b;
c2 = c2 * b;
while 1
    xnew = x+a*p;
    if f(xnew) <= fx + c1*a
        gxp = grad(xnew);
        if strong
            if abs(gxp'*p) <= abs(c2)
                break;
            end
        else
            if gxp'*p >= c2
                break;
            end
        end
    end
    a = rho * a;
end

function Goldstein(fun, x, p, a0, rho, c)
fx = fun(x);
a = a0;
fxp = f(x + a*p);
b = grad'*p;
c1 = b*c;
c2 = (1-c)*b;
while fxp > fx + c1*a || fxp < fx + c2*a
    a = rho * a;
    fxp = f(x + a*p);
end

function LineSearch()
