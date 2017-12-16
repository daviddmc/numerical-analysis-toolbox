function [A, J] = Jacobi(A, tol)
% Jacobi    Jacobi method for eigenvalue and eigenvector
%   E = Jacobi(A) produces a column vector E containing the eigenvalues 
%   of a symmetric (Hermitian) matrix A.
% 
%   [D, J] = Jacobi(A) produces a diagonal matrix D of eigenvalues and a 
%   unitary matrix J so that D = J'*A*J, columns of J are the corresponding
%   eigenvectors.
%
%   See also

%   Copyright 2017 Junshen Xu

if ~exist('tol','var')
    tol = 1e-10;
end

n = size(A, 1);
delta = Norm(A - diag(diag(A)), 'fro') / n;

if nargout>1
    J = eye(n);
end

while(delta > tol)
    if isreal(A)
        for k = 1 : n - 1
            for j = k+1 : n
                if(abs(A(j,k)) > delta)
                    tau = (A(j,j) - A(k,k)) / (2 * A(j,k));
                    t = sign(tau) / (abs(tau) + sqrt(1 + tau^2));
                    c = 1 / sqrt(1 + t^2);
                    s = c*t;
                    q = [c s; -s' c];
                    A([j, k], :) = q * A([j, k],:);
                    A(:, [j, k]) = A(:, [j, k]) * q';
                    if nargout > 1
                        J(:, [j, k]) = J(:, [j, k]) * q';
                    end
                end
            end
        end
    else
        for k = 1 : n - 1
            for j = k+1 : n
                if(abs(A(j,k)) > delta)
                    d = real(A(j, j) - A(k, k));
                    if d > 0
                        dMax = d + sqrt(4*abs(A(j,k))^2 + d^2);
                    else
                        dMax = d - sqrt(4*abs(A(j,k))^2 + d^2);
                    end
                    t = min(abs(2*A(j,k)/dMax), 1);
                    c = 1 / sqrt(1+t^2);
                    s = c*t*sign(A(k,j));
                    q = [c s'; -s c];
                    A([j, k], :) = q' * A([j, k],:);
                    A(:, [j, k]) = A(:, [j, k]) * q;
                    if nargout > 1
                        J(:, [j, k]) = J(:, [j, k]) * q;
                    end
                end
            end
        end
    end
    delta = delta / n;
end

if nargout < 2
    A = real(diag(A));
end




