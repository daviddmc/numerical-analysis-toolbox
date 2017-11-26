function A = Inverse( A )
%INVERSE Summary of this function goes here
%   Detailed explanation goes here

n = size(A, 1);
CheckSquareMatrix(A);

for k = 1 : n
    
    A(k, k) = 1 / A(k, k);
    for j = [1:k-1 k+1:n]
        A(k, j) = A(k, j) * A(k, k);
    end
    
    for i = [1:k-1 k+1:n]
        for j = [1:k-1 k+1:n]
            A(i,j) = A(i,j) - A(i,k) * A(k,j);
        end
    end
    
    for i = [1:k-1 k+1:n]
        A(i, k) = -A(i,k) * A(k, k);
    end
    
end

