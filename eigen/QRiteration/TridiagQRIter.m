function [D, Q] = TridiagQRIter(D, Q, tol)
% TridiagQRIter   QR iteration for symmetric tridiagonal marix.
%   D = TridiagQRIter(T) is the diagonal form of the symmetric tridiagonal
%   matrix T, whose elements are the eigenvalue of T. Wilkinson shift is
%   used in this algorithm.
%
%   [D, Q] = TridiagQRIter(T) produces a diagonal matrix D and a 
%   unitary matrix Q so that D = Q*T*Q', columns of Q are the 
%   corresponding eigenvectors for eigenvalues of T. 
% 
%   [D, Q] = TridiagQRIter(T, P) assumes that P is a unitary matrix of the
%   same size as T and produces a diagonal matrix D and a unitary matrix Q 
%   so that D = Q*P'*T*P*Q', columns of Q is the corresponding eigenvectors
%   for eigenvalues of P'*T*P.
%
%   [D, Q] = TridiagQRIter(T, P, TOL) specifies the tolerance of QR
%   iteration. If TOL is [] the default value, 1e-10, will be used.
%
%   See also SingleHessenbergQRIter, DoubleHessenbergQRIter, BidiagQRITer,
%   SymEigen.

%   Copyright 2017 Junshen Xu

flagQ = nargout > 1;
n = size(D, 1);
if flagQ && ~exist('Q', 'var')
    Q = eye(n);
end

if ~exist('tol', 'var')
    tol = 1e-10;
end

while(1)
    D(1:n+1:end) = real(D(1:n+1:end));
    for i = 1:n-1
        if abs(D(i+1,i)) < tol*(abs(D(i,i) + abs(D(i+1,i+1))))
            D(i+1,i) = 0;
            D(i,i+1) = 0;
        end
    end
    [i, j] = FindUnreduced(D, n);
    if j == 1
        break;
    end
    [x, z] = WilkinsonShift(D, i, j);
    if(z ~= 0)
        [c, s] = GivensRotation(x, z);
        m = min(i+2, j);
        D([i, i+1], i:m) = [c, s; -s', c] * D([i, i+1], i:m);
        D(i:m, [i, i+1]) = D(i:m, [i, i+1]) * [c, -s; s', c];
        if(flagQ)
            Q(:, [i, i+1]) = Q(:,[i,i+1])* [c, -s; s', c];
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
        D([k, k+1], k-1:m) = [c, s; -s', c] * D([k, k+1], k-1:m);
        D(k-1:m, [k, k+1]) = D(k-1:m, [k, k+1]) * [c, -s; s', c];
        D(k+1,k-1) = 0;
        D(k-1,k+1) = 0;
        if(flagQ)
            Q(:, [k, k+1]) = Q(:,[k,k+1])* [c, -s; s', c];
        end
    end
end

end

function [i, j] = FindUnreduced(D, n)
    for j = n:-1:1
        if j > 1
            if(D(j-1, j)~= 0)
                break;
            end
        end
    end
    for i = j-1:-1:1
        if i > 1
            if(D(i-1, i) == 0)
                break;
            end
        end
    end
end

function [x, z] = WilkinsonShift(D, i, j)
    d = (D(j-1,j-1) - D(j,j))/2;
    c = D(j,j-1)'*D(j,j-1);
    if d == 0
        mu = D(j,j) - sqrt(c);
    else
        mu = D(j,j) - c/(d+sign(d)*sqrt(d^2+c));
    end
    x = D(i,i) - mu;
    z = D(i+1,i);
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