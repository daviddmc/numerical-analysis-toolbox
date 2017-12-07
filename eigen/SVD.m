function [S, U, V] = SVD( A, econ )
% SVD    Singular value decomposition.
%   [S, U, V] = SVD(X) produces a diagonal matrix S, of the same dimension 
%   as X and with nonnegative diagonal elements, and unitary matrices U and 
%   V so that X = U*S*V'.
%
%   S = svd(X) returns a vector containing the singular values.
%   
%   [S, U, V] = SVD(X, 0) also produces the "economy size" decomposition. 
%   If X is m-by-n with m > n, then only the first n columns of U are 
%   computed and S is n-by-n. For m < n, only the first m columns of V are 
%   computed and S is m-by-m.
%
%   See also 

%   Copyright 2017 Junshen Xu

flagU = nargout > 1;
flagV = nargout > 2;

[m, n] = size(A);
if n > m
    flagT = 1;
    A = A';
else
    flagT = 0;
end

if(~exist('econ','var'))
    econ = 1;
end

if flagU
    if flagV
        [A, U, V] = BidiagReduction(A, econ); 
    else
        [A, U] = BidiagReduction(A, econ); 
    end
else
    A = BidiagReduction(A, econ); 
end
[m, n] = size(A);

epsilon = 1e-12;
Bnorm = max(abs(diag(A)) + [diag(A,1);0]);
tol = epsilon * Bnorm;

while 1
    for i = 1 : n - 1
        if(abs(A(i, i+1)) < epsilon * (abs(A(i, i)) + abs(A(i+1, i+1))))
            A(i, i+1) = 0;
        end
    end
    [i,j] = FindUnreduced(A, n);
    if(j == 1)
        break;
    end
    flagDiagZero = 0;
    for k = i : j-1
        if(abs(A(i, i)) < tol)
            A(i, i) = 0;
            A(i, i+1) = 0;
            flagDiagZero = 1;
        end
    end
    if flagDiagZero
        continue;
    end
    
    [y, z] = GolubKahan(A, i, j);
    for k = i : j-1
        if(z ~= 0)
            [c, s] = GivensRotation(y, z);
            q = max(k-1, 1);
            G = [c, s;-s' c];
            A(q:k+1,[k,k+1]) = A(q:k+1,[k, k+1]) * G;
            if k > i
                A(k-1,k+1) = 0;
            end
            if flagV
                V(:, [k,k+1]) = V(:, [k,k+1]) * G;
            end
            y = A(k, k);
            z = A(k+1, k);
            if(z~=0)
                [c, s] = GivensRotation(y, z);
                q = min(k+2, n);
                GT = [c -s';s c];
                A([k, k+1], k:q) = GT * A([k,k+1], k:q);
                A(k+1, k) = 0;
                if flagU
                    U(:, [k,k+1]) = U(:, [k,k+1]) * GT';
                end
                if k < j - 1
                    y = A(k, k+1);
                    z = A(k ,k+2);
                end
            end
        end
    end 
end


if flagV
    V = V * diag(sign(diag(A)'));
end

if ~flagU
    S = abs(diag(A));
else
    S = zeros(size(A));
    S(1:m+1:end) = abs(diag(A));
end

if flagT
    S = S';
    tmp = V;
    V = U;
    U = tmp;
end

end

function [i,j] = FindUnreduced(A, n)
    for j = n : -1 : 1
        if(j > 1)
            if(A(j-1,j)~=0)
                break;
            end
        end
    end
    for i = j-1 : -1 : 1
        if(i > 1)
            if(A(i-1,i) == 0)
                break;
            end
        end
    end
end

function [y, z] = GolubKahan(A, i, j)
    c = abs(A(j-1, j))^2 + abs(A(j, j))^2;
    if(j - i> 1)
        a = abs(A(j-2, j-1))^2 + abs(A(j-1, j-1))^2;
    else
        a = abs(A(j-1, j-1))^2;
    end
    b = abs(A(j-1, j)*A(j-1, j-1));
    if(a > c)
        mu = (a+ c - sqrt((a-c)^2+4*b^2)) / 2;
    else
        mu = (a+ c + sqrt((a-c)^2+4*b^2)) / 2;
    end
    y = A(i, i)^2 - mu;
    z = A(i, i)*A(i, i+1);
end

function [c, s] = GivensRotation(y, z)
    
    t = sqrt(abs(y)^2 + abs(z)^2);
    c = abs(y) / t;
    if c == 0
        s = -z / t;
    else
        s = -z * sign(y') / t;
    end
    
end

