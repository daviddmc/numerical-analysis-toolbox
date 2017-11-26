function [L, D] = Cholesky(A)
% Cholesky decomposition
% L * L^T = A or L * D * L^T = A
% input
% A
% output
% L
% D

CheckSquareMatrix(A, 'A');
n = size(A, 1);

if nargout > 1
    L = eye(n);
    D = zeros(n, 1);
    D(1) = A(1,1);
    T = zeros(size(A));
    for ii = 2 : n
        for jj = 1 : ii - 1
            T(ii, jj) = A(ii, jj) - T(ii, 1:jj-1) * L(jj, 1:jj-1)';
            L(ii, jj) = T(ii, jj) / D(jj);
        end
        D(ii) = A(ii, ii) - T(ii, 1:ii-1) * L(ii, 1:ii-1)';
    end
else
    L = zeros(size(A));
    for jj = 1 : n
        L(jj, jj) = sqrt(A(jj, jj) - sum(L(jj,1:jj-1).^2));
        for ii = jj+1 : n
            L(ii, jj) = (A(ii, jj) - L(ii,1:jj-1) * L(jj,1:jj-1)') / L(jj,jj);
        end
    end
end


