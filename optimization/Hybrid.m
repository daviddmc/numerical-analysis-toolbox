function [ x, flag, iter ] = Hybrid( fun, jac, x0, maxIter, tol1, tol2, tau )
%HYBRID Summary of this function goes here
%   Detailed explanation goes here

n = length(x0);
B = eye(n);
I = eye(n);
x = x0;
%mu = mu0;
nu = 2;
f = fun(x);
F = f'*f / 2;
J = jac(x);
A = J'*J;
mu = tau * max(diag(A));
g = J'*f;
gInfNorm = max(abs(g));
found = gInfNorm < tol1;
flag = 1;
count = 0;
method = 1; % L-M
e = sqrt(eps);
for iter = 1 : maxIter
    if found
        flag = 0;
        break
    end
    if(method) % L-M
        xNew = x;
        %h = (A+mu*I)\(-g);
        h = SPDSolve(A+mu*I, -g);
        if Norm(h) < tol2 * (Norm(x) + tol2)
            found = 1;
        else
            xNew = x + h;
            fNew = fun(xNew);
            FNew = fNew'*fNew / 2;
            JNew = jac(xNew);
            gNew = JNew'*fNew;
            rho = (F - FNew) / (h'*(mu*h-g)/2);
            if rho > 0
                better = 1;
                
                gNewInfNorm = max(abs(gNew));
                found = gNewInfNorm < tol1;
                if gNewInfNorm < 0.02 * FNew
                    count = count + 1;
                    if count == 3
                        method = 0; % Q-N
                        % switch
                        delta = max([1.5*tol2*(Norm(x) + tol2),...
                            Norm(h) / 5]);
                        count = 0;
                    end
                else
                    count = 0;
                end
                mu = mu * max([1/3, 1-(2*rho-1)^3]);
                nu = 2;
            else
                mu = mu * nu;
                nu = 2 * nu;
                count = 0;
                better = 0;
            end       
        end
    else % Quasi-Newton
        xNew = x;
        better = 0;
        %h = B\(-g);
        h = SPDSolve(B, -g);
        hNorm = Norm(h);
        if hNorm < tol2 * (Norm(x) + tol2)
            found = 1;
        else
            if hNorm > delta
                h = delta / hNorm * h;
            end
            
            xNew = x + h;
            fNew = fun(xNew);
            FNew = fNew'*fNew / 2;
            JNew = jac(xNew);
            gNew = JNew'*fNew;
            gNewInfNorm = max(abs(gNew));
            if gNewInfNorm < tol1
                found = 1;
            else
                better = (FNew < F) || ...
                    (FNew < (1+e)*F && gNewInfNorm < gInfNorm);
                if gNewInfNorm >= gInfNorm
                    method = 1;
                else
                    % update delta
                    L = -h'*g - h'*A*h / 2;
                    rho = (F - FNew) / L;
                    if rho < 0.25
                        delta = delta / 2;
                    elseif rho > 0.75
                        delta = max([delta, 3*hNorm]);
                    end
                end
            end
        end
    end
    
    
    if ~found
        %update B
        ANew = JNew'*JNew;
        y = ANew*h + gNew - J'*fNew;
        hy = h'*y;
        if hy > 0
            v = B*h;
            B = B + (y/hy)*y' - (v/(h'*v))*v';
        end
        
        if better
            x = xNew;
            g = gNew;
            gInfNorm = gNewInfNorm;
            F = FNew;
            J = JNew;
            A = ANew;
        end
    end
    %
    
end

