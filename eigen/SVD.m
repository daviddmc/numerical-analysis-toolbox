function S = SVD( A )
%SVD Summary of this function goes here
%   Detailed explanation goes here

A = BidiagonalReduction( A ); 
Bd = diag(A);
Bs = diag(A, 1);
n = length(Bd);

epsilon = 1e-10;
Bnorm = max(abs(Bd) + [abs(Bs);0]);
tol = epsilon * Bnorm;

while 1
    for ii = 1 : n - 1
        if(abs(Bs(ii)) < epsilon * (abs(Bd(ii)) + abs(Bd(ii+1))))
            Bs(ii) = 0;
        end
    end
    for ii = 1 : n
        if(abs(Bd(ii)) < tol)
            Bd(ii) = 0;
            if(ii < n)
                Bs(ii) = 0;
            end
        end
    end
    for jj = n : -1 : 1
        if(jj > 1)
            if(Bs(jj-1)~=0)
                break;
            end
        end
    end
    if(jj == 1)
        break;
    end
    
    for ii = jj-1 : -1 : 1
        if(ii > 1)
            if(Bs(ii-1) == 0)
                break;
            end
        end
    end
    
    c = Bs(jj-1)^2 + Bd(jj)^2;
    if(jj - ii> 1)
        a = Bs(jj-2)^2 + Bd(jj-1)^2;
    else
        a = Bd(jj-1)^2;
    end
    b = Bs(jj-1)*Bd(jj-1);
    if(a > c)
        mu = (a+ c - sqrt((a-c)^2+4*b^2)) / 2;
    else
        mu = (a+ c + sqrt((a-c)^2+4*b^2)) / 2;
    end
    y = Bd(ii)^2 - mu;
    z = Bd(ii)*Bs(ii);
    for k = ii : jj-1
        if(z ~= 0)
            t = sqrt(abs(y)^2 + abs(z)^2);
            c = y/t;
            s = -z/t;
            tmp = Bd(k);
            Bd(k) = c * Bd(k) - s * Bs(k);
            Bs(k) = s * tmp + c * Bs(k);
            if(k > ii)
                Bs(k-1) = Bs(k-1) * c - s * z;
            end
            z = - Bd(k+1) * s;
            y = Bd(k);
            Bd(k+1) = c * Bd(k+1);

            if(z~=0)
                t = sqrt(abs(y)^2 + abs(z)^2);
                c = y/t;
                s = -z/t;
                Bd(k) = c * Bd(k) - s * z;
                tmp = Bs(k);
                Bs(k) = c * Bs(k) - s * Bd(k+1);
                Bd(k+1) = s * tmp + c * Bd(k+1);
                if k < jj - 1
                    y = Bs(k);
                    z = -s * Bs(k+1);
                    Bs(k+1) = c * Bs(k+1);
                end
            end
        end
    end 
end
        
S = abs(Bd);

