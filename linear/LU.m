function [L, U, P] = LU(A)
% LU decomposition with partial pivot
% PA = LU
% input
% A
% output
% L
% U
% P

CheckSquareMatrix(A, 'A');
n = size(A, 1);
L = eye(n);
U = A;
p = zeros(2 ,n - 1);
pp = 1:n;
P = eye(n);

for c = 1 : n-1
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

for c = 1 : n-1
    for k = c+1:n-1
        L([p(1, k), p(2, k)], c) = L([p(2,k), p(1,k)], c);
    end
end

P = P(pp, :);
