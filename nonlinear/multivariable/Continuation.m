function x = Continuation( fun, jac, N, x0)
%BROYDEN Summary of this function goes here
%   Detailed explanation goes here

h = 1/N;
b = -h*fun(x0);
x = x0;

for k = 1:N
    k1 = GaussianElimination(jac(x), b);
    k2 = GaussianElimination(jac(x + 0.5*k1), b);
    k3 = GaussianElimination(jac(x + 0.5*k2), b);
    k4 = GaussianElimination(jac(x + k3), b);
    x = x + (k1+2*k2+2*k3+k4) / 6;
end

end

