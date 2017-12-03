function [L, D] = LDL( A )
% LDL   ldl' decomposition for Hermitian matrix.
%   Assuming that A is Hermitian matrix (that is, A == A'),
%   [L, D] = LDL(A) stores a vector D and 
%   
%
%   See also

%   Copyright 2017 Junshen Xu

CheckSquareMatrix(A, 'A');
n = size(A, 1);

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

