function [L, D] = Cholesky(A, shape)
%Cholesky   Cholesky or LDL decomposition of positive definite matrix.
%   Cholesky(A) uses only the diagnoal and lower triangle of A. The upper 
%   triangle is assumed to be the (complex conjugate) transpose of the 
%   lower triangle. If A is positive definite, then L = Cholesky(A)
%   produces a lower triangular L so that L*L' = A. If A is not positive
%   definite, an error message is printed.
%
%   R = Cholesky(A, 'upper') uses ony the diagonal and the upper triangle
%   of A to produce a upper triangular R so that R'*R = A. If A is not
%   positive definite, an error message is printed.
%
%   [L, D] = Cholesky(A) produces a lower triangular L whose elements in 
%   main diagonal equal 1 and a positive vector D so that L*diag(D)*L' = A.
%   If A is not positive definite, an error message is printed. 
% 
%   [R, D] = Cholesky(A, 'upper') produces a upper triangular R whose 
%   elements in main diagonal equal 1 and a positive vector D so that 
%   R'*diag(D)*R = A. If A is not positive definite, an error message is 
%   printed. 
%
%   See also LU, 

%   Copyright 2017 Junshen Xu

CheckSquareMatrix(A, 'A');
n = size(A, 1);

D = diag(A);
if any(~(abs(imag(D)) < eps) | real(D) < eps)
    error('Matrix must be positive definite.');
end
    
if(~exist('shape','var'))
    shape = 'lower';
end

if nargout > 1
    diagnoalIdx = 1:n+1:n^2;
    if(strcmpi(shape, 'lower'))
        for j=1:n
            v = A(j,1:j-1) .*  A(diagnoalIdx(1:j-1));
            v = v';
            A(j,j) = A(j,j) - real(A(j,1:j-1)*v);
            if(A(j,j) < eps)
                error('Matrix must be positive definite.');
            end
            A(j+1:n,j) = (A(j+1:n,j) - A(j+1:n,1:j-1)*v) / A(j,j);
        end
        L = tril(A, -1);
        L(diagnoalIdx) = 1;
        D = diag(A);
    elseif(strcmpi(shape, 'upper'))
        for j=1:n
            v = A(1:j-1,j) .*  A(diagnoalIdx(1:j-1)');
            v = v';
            A(j,j) = A(j,j) - real(v*A(1:j-1,j));
            if(A(j,j) < eps)
                error('Matrix must be positive definite.');
            end
            A(j,j+1:n) = (A(j,j+1:n) - v*A(1:j-1,j+1:n)) / A(j,j);
        end
        L = triu(A, 1);
        L(diagnoalIdx) = 1;
        D = diag(A);
    else
        error('Shape flag must be ''upper'' or ''lower''.');
    end
    return
end

if(strcmpi(shape, 'lower'))
    for jj = 1 : n
        if jj > 1
            A(jj:n, jj) = A(jj:n, jj) - A(jj:n, 1:jj-1) * A(jj, 1:jj-1)';
        end
        A(jj, jj) = real(A(jj,jj));
        if(A(jj, jj) <= eps)
            error('Matrix must be positive definite.');
        end
        A(jj:n,jj) = A(jj:n, jj) / sqrt(A(jj,jj));
        L = tril(A);
    end
elseif(strcmpi(shape, 'upper'))
    for jj = 1 : n
        if jj > 1
            A(jj,jj:n) = A(jj,jj:n) - A(1:jj-1,jj)' * A(1:jj-1,jj:n);
        end
        A(jj, jj) = real(A(jj,jj));
        if(A(jj, jj) <= eps)
            error('Matrix must be positive definite.');
        end
        A(jj,jj:n) = A(jj,jj:n) / sqrt(A(jj,jj));
        L = triu(A);
    end
else
    error('Shape flag must be ''upper'' or ''lower''.');
end


