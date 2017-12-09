function x = Broyden( fun, x0 , tol)
%BROYDEN Summary of this function goes here
%   Detailed explanation goes here

n = length(x0);
fx = fun(x0);
Binv = eye(n);
x = x0;
while Norm(fx) > tol
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

