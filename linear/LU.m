function [L, U, P, Q] = LU(A, pType)
% LU    lu decomposition.
%   [L,U] = LU(A) stores an upper triangular matrix in U and a
%   "psychologically lower triangular matrix" (i.e. a product of lower
%   triangular and permutation matrices) in L, so that A = L*U. A can be
%   rectangular.
%
%   [L, U, P] = LU(A) returns unit lower triangular matrix L, upper
%   triangular matrix U, and permutation matrix P so that P*A = L*U.
%
%   [L, U, p] = LU(A, 'vector') returns the permutation information as a
%   vector instead of a matrix.  That is, p is a row vector such that
%   A(p,:) = L*U.  Similarly, [L, U, P] = LU(A, 'matrix') returns a
%   permutation matrix P.  This is the default behavior.
%
%   [L, U, P, Q] = LU(A) returns unit lower triangular matrix L, upper
%   triangular matrix U, and two permutation matrix P and Q so that 
%   P*A*Q = L*U. The column permutation Q is chosen so that ABS(DIAG(U)) 
%   is decreasing.
%
%   [L, U, p, q] = LU(A,'vector') returns two row vectors p and q so that
%   A(p,q) = L*U.  Using 'matrix' in place of 'vector' returns permutation
%   matrices.
%
%   See also Cholesky, QR.

%   Copyright 2017 Junshen Xu

if(~ismatrix(A) || ~isnumeric(A))
    error('A should be a numeric matrix')
end

if(~exist('pType', 'var'))
    pType = 'matrix';
end

[m, n]= size(A);
K = min(m, n);

L = eye(m, n);
U = A;
p = 1:m;

if nargout > 3
    q = 1:n;
    for c = 1 : K
        % find pivot
        [M, idx1] = max(abs(U(c:end, c:end)));
        [~, idx2] = max(M);
        idx1 = idx1(idx2) + c - 1;
        idx2 = idx2 + c - 1;
        U([idx1 ,c], c:end) = U([c, idx1], c:end);
        U(:, [idx2, c]) = U(:, [c, idx2]);
        % update p
        p([idx1, c]) = p([c, idx1]);
        q([idx2, c]) = q([c, idx2]);
        % update L
        if(abs(U(c,c)) < eps)
            break;
        end
        l = U(c+1:end,c)/ U(c,c);
        L(c+1:end, c) = l;
        L([idx1, c], 1:c-1) = L([c, idx1], 1:c-1);
        % update U
        U(c+1:end, c:end) = U(c+1:end, c:end) - bsxfun(@times, l, U(c, c:end)) ;
    end
else
    for c = 1 : K
        % find pivot
        [~, idx] = max(abs(U(c:end, c)));
        idx = idx + c - 1;
        U([idx ,c], c:end) = U([c, idx], c:end);
        % update p
        p([idx, c]) = p([c, idx]);
        % update L
        if(abs(U(c,c)) < eps)
            continue;
        end
        l = U(c+1:end,c)/ U(c,c);
        L(c+1:end, c) = l;
        L([idx, c], 1:c-1) = L([c, idx], 1:c-1);
        % update U
        U(c+1:end, c:end) = U(c+1:end, c:end) - bsxfun(@times, l, U(c, c:end)) ;
    end
end

if(m > n)
    U = U(1:n, :);
else
    L = L(:, 1:m);
end

if(nargout > 2)
    if(strncmpi(pType, 'v', 1))
        P = p;
    else
        P = eye(m);
        P = P(p, :);
    end
    if(nargout > 3)
        if(strncmpi(pType, 'v',1))
            Q = q;
        else
            Q = eye(n);
            Q = Q(:, q);
        end
    end
else
    L(p, :) = L;
end
