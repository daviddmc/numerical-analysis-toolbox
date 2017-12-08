function x = SR1TrustRegion(fun, grad, delta0, tol, eta, r,x0 )
%SR1TRUSTREGION Summary of this function goes here
%   Detailed explanation goes here

eta = 1e-3;
r = 1e-8;
delta = delta0;
x = x0;
g = grad(x0);
f = fun(x0);
B = eye(length(x0));
while Norm(g) > tol
    %s = DogLeg(g, delta, B);
    s = CauchyPoint( g, delta, B );
    sNorm = Norm(s);
    xNew = x + s;
    fNew = fun(xNew);
    gNew = grad(xNew);
    y = gNew - g;
    ared = f - fNew;
    Bs = B*s;
    pred = - (g'*s + s'*Bs / 2);
    ratio = ared / pred;
    if ratio > eta
        x = xNew;
        g = gNew;
        f = fNew;
    end
    if ratio > 0.75
        if sNorm > 0.8*delta
            delta = 2*delta;
        end
    elseif ratio < 0.1
        delta = delta / 2;
    end
    yBs = y - Bs;
    syBs = s' * yBs;
    if abs(syBs) >= r*sNorm*Norm(yBs)
        B = B + (yBs)*(yBs)' / syBs;
    end
end

