function A = Inverse( A )
% Inverse   Matrix inverse.
% 
%   Inverse(A) is the inverse of the square matrix A calcualted by
%   Gauss-Jordan method with complete pivoting.
%
%   See also

%   Copyright 2017 Junshen Xu

CheckSquareMatrix(A);
n = size(A, 1);
rowSwap = zeros(2, n);
colSwap = zeros(2, n);

for k = 1 : n
    
    [m, idx1] = max(A(k:end, k:end));
    [~, idx2] = max(m);
    rowSwap(:, k) = [idx1(idx2) + k - 1, k];
    colSwap(:, k) = [idx2 + k - 1, k];
    
    A([rowSwap(1,k) rowSwap(2,k)], :) = A([rowSwap(2,k) rowSwap(1,k)], :);
    A(:, [colSwap(1,k) colSwap(2,k)]) = A(:, [colSwap(2,k) colSwap(1,k)]);
    
    if(A(k, k) == 0)
        error('Matrix is singular to working precision');
    end
    
    a = 1 / A(k, k);
    aC = A(k,:) * a;
    aR = A(:,k) * -a;
    A = A - A(:, k) * aC;
    A(k,:) = aC;
    A(:,k) = aR;
    A(k, k)= a;
    
end

for k = n:-1:1
    A(:, [rowSwap(1,k) rowSwap(2,k)]) = A(:, [rowSwap(2,k) rowSwap(1,k)]);
    A([colSwap(1,k) colSwap(2,k)], :) = A([colSwap(2,k) colSwap(1,k)], :);
end

