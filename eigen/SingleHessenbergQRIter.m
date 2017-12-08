function [D, Q] = SingleHessenbergQRIter(D, Q, tol)
%QRITER Summary of this function goes here
%   Detailed explanation goes here

flagQ = nargout > 1;
n = size(D, 1);
if ~exist('tol','var')
    tol = 1e-10;
end

while(1)
    for i = 1:n-1
        if abs(D(i+1,i)) < tol*(abs(D(i,i) + abs(D(i+1,i+1))))
            D(i+1,i) = 0;
        end
    end
    [i, j] = FindUnreduced(D, n);
    if j == 1
        break;
    end
    [x, z] = SingleShift(D, i, j);
    if(z ~= 0)
        [c, s] = GivensRotation(x, z);
        m = min(i+2, j);
        if flagQ
            D([i, i+1], i:end) = [c, s; -s', c] * D([i, i+1], i:end);
            D(1:m, [i, i+1]) = D(1:m, [i, i+1]) * [c, -s; s', c];
            Q(:, [i, i+1]) = Q(:,[i,i+1])* [c, -s; s', c];
        else
            D([i, i+1], i:j) = [c, s; -s', c] * D([i, i+1], i:j);
            D(i:m, [i, i+1]) = D(i:m, [i, i+1]) * [c, -s; s', c];
        end
    end
    
    for k = i+1:j-1
        x = D(k, k-1);
        z = D(k+1, k-1);
        if(z == 0)
            continue;
        end
        [c, s] = GivensRotation(x, z);
        m = min(k+2, j);
        if flagQ
            D([k, k+1], k-1:end) = [c, s; -s', c] * D([k, k+1], k-1:end);
            D(1:m, [k, k+1]) = D(1:m, [k, k+1]) * [c, -s; s', c];
            Q(:, [k, k+1]) = Q(:,[k,k+1])* [c, -s; s', c];
        else
            D([k, k+1], k-1:j) = [c, s; -s', c] * D([k, k+1], k-1:j);
            D(i:m, [k, k+1]) = D(i:m, [k, k+1]) * [c, -s; s', c];
        end
    end
end

end

function [i, j] = FindUnreduced(D, n)
    for j = n:-1:1
        if j > 1
            if(D(j, j-1)~= 0)
                break;
            end
        end
    end
    for i = j-1:-1:1
        if i > 1
            if(D(i, i-1) == 0)
                break;
            end
        end
    end
end

function [x, z] = SingleShift(D, i, j)
    if j - i > 1 || imag(D(j,j)) ~= 0
        x = D(i, i) - D(j, j);
    else
        ad = D(j-1, j-1) - D(j,j);
        bc = D(j,j-1)*D(j-1,j);
        mu = (ad + sqrt(ad^2 + 4*bc)) / 2;
        x = D(i,i) - mu;
    end
    z = D(i+1, i);
end

function [c, s] = GivensRotation(x, z)
    t = sqrt(abs(x)^2 + abs(z)^2);
    c = abs(x) / t;
    if c == 0
        s = z' / t;
    else
        s = z' * sign(x) / t;
    end
end