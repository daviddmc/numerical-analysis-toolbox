function [x, flag, iter, delta] = LSTrustRegionDogLeg( fun, jac, x0, maxIter, delta0,tol1, tol2, tol3 )
%LSTRUSTREGION Summary of this function goes here
%   Detailed explanation goes here

x = x0;
delta = delta0;
f = fun(x);
F = f'*f / 2;
J = jac(x);
[n,m] = size(J);
g = J'*f;
flag = 1;
found = max(abs(f)) < tol3 || max(abs(g)) < tol1;
for iter = 1 : maxIter
    if found
        flag = 0;
        break
    end
    Jg = J*g;
    gg = g'*g;
    alpha = gg / (Jg'*Jg);
    hsd = -alpha*g;
    hsdNorm = Norm(hsd);
    if n == m
        hgn = GaussianElimination(J, -f);
    else
        hgn = LS(J, -f);
    end
    
    if Norm(hgn) <= delta
        hdl = hgn;
        L = F;
    elseif hsdNorm >= delta
        hdl = delta / hsdNorm * hsd;
        L = delta * (2*hsdNorm -delta) / (2*alpha);
    else
        a2 = hsdNorm^2;
        ba = hgn - hsd;
        ba2 = ba'*ba;
        c = hsd' * ba;
        if c <= 0
            beta = (-c + sqrt(c^2 + ba2*(delta^2-a2))) / ba2;
        else
            beta = (delta^2 - a2) / (c + sqrt(c^2 + ba2*(delta^2 - a2)));
        end
        hdl = hsd + beta*ba;
        L = 0.5*alpha*(1-beta)^2*gg + beta*(2-beta)*F;
    end
    
    if Norm(hdl) < tol2*(Norm(x) + tol2)
        found = 1;
    else
        xNew = x + hdl;
        fNew = fun(xNew);
        FNew = fNew'*fNew / 2;
        rho = (F - FNew) / L;
        if rho > 0
            x = xNew;
            f = fNew;
            F = FNew;
            J = jac(x);
            g = J'*f;
            found = max(abs(f)) < tol3 || max(abs(g)) < tol1;
        end
        
        if rho > 0.75
            delta = max([delta, 3*Norm(hdl)]);
        elseif rho < 0.25
            delta = delta / 2;
            found = delta < tol2*(Norm(x) + tol2);
        end
    end
end


end

