function A = Inverse(A, method)
% Inverse   Matrix inverse.
% 
%   B = Inverse(A) is the inverse of the square matrix A.
%
%   B = Inverse(A, method) specifies the algorithm for matrix inversion. 
%   The available method are:
%
%     'GaussJordan' - (default) Gauss-Jordan elimiation
%     'LU'          - LU decomposition
%     'QR'          - QR decomposition
%     'SVD'         - singular value decomposition    
%
%   See also

%   Copyright 2017 Junshen Xu

CheckSquareMatrix(A);
n = size(A, 1);

if(~exist('method','var') || isemtpy(method))
    method = 'gauss';
end

if strncmpi(method, 'lu', 1)

    [L, U, p] = LU(A, 'vector');
    A = TriangleInverse(U, 'upper') * TriangleInverse(L, 'lower');
    A(:, p) = A;
    
elseif strncmpi(method, 'qr', 1)
    [Q, R] = QR(A);
    A = TriangleInverse(R, 'upper') * Q';

elseif strncmpi(method, 'svd', 1)
    error('not implemented yet');
elseif strncmpi(method, 'GaussJordan', 1)

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
else
    error('method error');
end

