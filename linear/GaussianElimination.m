function X = GaussianElimination(A, B, pivot)
% GaussianElimination Solve linear system using Gaussian elimination.
%   X = GaussianElimination(A, B) solves the linear system A*X=B using
%   Gaussian elimination with partial pivoting.
%
%   X = GaussianElimination(A, B, PIVOT) specifies the pivoting method.
%   The available methods are:
%   
%     'none'     - use no pivot
%     'partial'  - (default) use the maximum in the column as pivot
%     'complete' - use the maximnu in the sub-matrix as pivot
%
%   See also

%   Copyright 2017 Junshen Xu

    CheckSquareMatrix(A);
    CheckMultiplicationSize(A,[],B);
    n = size(A, 1);
    
    if(~exist('pivot', 'var'))
        pivot = 'partial';
    end
    
    if(strcmpi(pivot, 'none'))
        pivotFlag = 0;
    elseif(strcmp(pivot, 'partial'))
        pivotFlag = 1;
    elseif(strcmp(pivot, 'complete'))
        pivotFlag = 2;
        columnPermute = 1 : n;
    else
        error('Invalid pivoting method.');
    end
        
    
    X = zeros(size(B));
    AB = [A B];
    % elimination
    for c = 1 : n-1
        if(pivotFlag == 1) % partial pivoting
            % find maximum
            [~, idx] = max(abs(AB(c:end,c)));
            idx = idx + c - 1;
            % swap
            AB([idx, c], c:end) = AB([c, idx], c:end);
        elseif(pivotFlag == 2) % complete pivoting
            % find maximum
            [m, idx1] = max(abs(AB(c:end,c:n)));
            [~, idx2] = max(m);
            idx1 = idx1(idx2) + c - 1;
            idx2 = idx2 + c - 1;
            % swap
            AB([idx1, c], c:end) = AB([c, idx1], c:end);
            AB(:, [idx2, c]) = AB(:, [c, idx2]);
            columnPermute([idx2, c]) = columnPermute([c, idx2]);
        end
        if(AB(c,c) == 0)
            if(pivotFlag)
                error('Matrix is singular to working precision.');
            else
                error('0 pivot found. Please try to use pivoting');
            end
        end
        AB(c+1:end, c:end) = AB(c+1:end, c:end) - bsxfun(@times, AB(c+1:end,c)/AB(c,c), AB(c, c:end));
    end
    % substitution
    X(n, :) = AB(n, n+1:end) / AB(n,n);
    for i = n-1 : -1 : 1
        X(i, :) = (AB(i, n+1:end) - AB(i,i+1:n) * X(i+1:n, :)) / AB(i,i);
    end
    
    if(pivotFlag == 2) 
        % If complete pivoting is used the order of x may be changed.
        X(columnPermute, :) = X;
    end
   
end