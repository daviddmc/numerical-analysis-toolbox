function x = GaussianElimination(A, b)
% Gaussian elimination with partial pivoting
% solution linear system Ax = b,
% where A is a square matrix and b is a vector.
% input
% A
% b
% output
% x
    
    CheckSquareMatrix(A);
    if(~isvector(b))
        error('b should be a vector');
    end
    if(~iscolumn(b))
        b = b(:);
    end
    CheckMultiplicationSize(A, b);

    n = length(b);
    x = zeros(n,1);
    A_b = [A b];
    % elimination
    for c = 1 : n-1
        A_b = FindPivot(A_b, c);
        A_b(c+1:end, c:end) = A_b(c+1:end, c:end) - bsxfun(@times, A_b(c+1:end,c), A_b(c, c:end)) / A_b(c,c);
    end
    % substitution
    x(n) = A_b(n, n+1)/A_b(n,n);
    for i = n-1 : -1 : 1
        x(i) = (A_b(i, n+1) - A_b(i,i+1:n) * x(i+1:n))/A_b(i,i);
    end
   
end

function x = FindPivot(x, c)
    % find maximum
    [~, idx] = max(abs(x(c:end,c)));
    idx = idx + c - 1;
    % swap
    a = x(idx, :);
    x(idx, :) = x(c, :);
    x(c, :) = a;
end