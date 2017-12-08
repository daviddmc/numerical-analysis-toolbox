function x = BFGS(fun, grad, tol, x0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

x = x0;
g = grad(x0);
H = eye(length(x0));
iter = 0;
while Norm(g) > tol
    iter = iter + 1;
    p = -(H*g);
    a = LineSearch(fun, grad, x, p, 1e-4, 0.9, 1, 10);
    s = a*p;
    xnew = x + s;
    gnew = grad(xnew);
    y = gnew - g;
    if iter == 1
        H = (y'*s)/(y'*y)*H;
    end
    rho = 1 / (y'*s);
    H = H - (rho*s)*(y'*H);
    H = H - (rho*(H*y))*s';
    H = H + rho*s*s';
    g = gnew;
    x = xnew;
end

