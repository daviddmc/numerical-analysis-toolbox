function [A, U, V] = BidiagQRIter( A, U, V, tol)
% BidiagQRIter   QR iteration for bidiagonal marix.
%   S = TBidiagQRIter(B) is the diagonal form of the bidiagonal matrix
%   matrix B of size M-by-N and M >= N. Golub - Kahan method is used in 
%   this algorithm.
%
%   [S, U] = BidiagQRIter(B, U0) assumes U0 is a matrix with orthogonal 
%   columns and produces a diagonal matrix S and a matrix U with orthogonal
%   columns.
%
%   [S, U, V] = BidiagQRIter(B, U0, V0) assumes U0 adn V0 are matrices with
%   orthogonal columns and produces a diagonal matrix S and matrices U and 
%   V with orthogonal columns.
%
%   [S, U, V] = BidiagQRIter(B, U0, V0, TOL) specifies the tolerance of QR
%   iteration. If TOL is [] the default value, 1e-10, will be used.
%
%   See also SingleHessenbergQRIter, DoubleHessenbergQRIter, 
%   TriidiagQRITer, SVD

%   Copyright 2017 Junshen Xu

flagU = nargout > 1;
flagV = nargout > 2;

if nargin < nargout
    error(' ');
end

n = size(A, 2);

if ~exist('tol','var')
    tol = 1e-10;
end

Bnorm = max(abs(diag(A)) + [diag(A,1);0]);
epsilon = tol * Bnorm;

while 1
    for i = 1 : n - 1
        if(abs(A(i, i+1)) < tol * (abs(A(i, i)) + abs(A(i+1, i+1))))
            A(i, i+1) = 0;
        end
    end
    [i,j] = FindUnreduced(A, n);
    if(j == 1)
        break;
    end
    flagDiagZero = 0;
    for k = i : j-1
        if(abs(A(i, i)) < epsilon)
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

