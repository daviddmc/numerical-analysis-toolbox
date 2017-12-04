function [L, U, P] = LU(A)
% LU    lu decomposition.
%   [L,U] = lu(A) stores an upper triangular matrix in U and a
%   "psychologically lower triangular matrix" (i.e. a product of lower
%   triangular and permutation matrices) in L, so that A = L*U. A can be
%   rectangular.
%
%   [L, U, P] = LU(A) returns unit lower triangular matrix L, upper
%   triangular matrix U, and permutation matrix P so that P*A = L*U.
%
%   See also Cholesky

%   Copyright 2017 Junshen Xu

if(~ismatrix(A) || ~isnumeric(A))
    error('A should be a numeric matrix')
end
[m, n]= size(A);
K = min(m, n);

L = eye(m, n);
U = A;
p = zeros(2 ,K);
pp = 1:m;
P = eye(m);

for c = 1 : K
    % find pivot
    [~, idx] = max(abs(U(c:end, c)));
    idx = idx + c - 1;
    U([idx ,c], c:end) = U([c, idx], c:end);
    % update p
    p(:, c) = [idx, c];
    pp([idx, c]) = pp([c, idx]);
    % update L
    l = U(c+1:end,c)/ U(c,c);
    L(c+1:end, c) = l;
    % update U
    U(c+1:end, c:end) = U(c+1:end, c:end) - bsxfun(@times, l, U(c, c:end)) ;
end

for c = 1 : K
    for k = c+1: K
        L([p(1, k), p(2, k)], c) = L([p(2,k), p(1,k)], c);
    end
end

if(m > n)
    U = U(1:n, :);
else
    L = L(:, 1:m);
end

P = P(pp, :);

if nargout < 3
    L = P' * L;
end
