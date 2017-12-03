function L = Cholesky(A, shape)
%Cholesky   Cholesky decomposition.
%   Cholesky(A) uses only the diagnoal and lower triangle of A. The upper 
%   triangle is assumed to be the (complex conjugate) transpose of the 
%   lower triangle. If A is positive definite, then L = Cholesky(A)
%   produces and lower triangular L so that L*L' = A. If A is not positive
%   definite, an error message is printed.
%
%   R = Cholesky(A, 'upper') uses ony the diagonal and the upper triangle
%   of A to produce a upper triangular R so that R'*R = A. If A is not
%   positive definite, an error message is printed.
% 
%   See also

%   Copyright 2017 Junshen Xu

CheckSquareMatrix(A, 'A');
n = size(A, 1);

if(~exist('shape','var'))
    shape = 'lower';
end

L = zeros(size(A));

if(strcmpi(shape, 'lower'))
    for jj = 1 : n
        temp = A(jj, jj) - sum(abs(L(jj,1:jj-1)).^2);
        if(~isreal(temp) || temp <= 0)
            error('Matrix must be positive definite.')
        else
            L(jj, jj) = sqrt(temp);
        end
        for ii = jj+1 : n
            L(ii, jj) = (A(ii, jj) - L(ii,1:jj-1) * L(jj,1:jj-1)') / L(jj,jj);
        end
    end
elseif(strcmpi(shape, 'upper'))
    for jj = 1 : n
        temp = A(jj, jj) - sum(abs(L(1:jj-1,jj)).^2);
        if(~isreal(temp) || temp <= 0)
            error('Matrix must be positive definite.')
        else
            L(jj, jj) = sqrt(temp);
        end
        for ii = jj+1 : n
            L(jj, ii) = (A(jj, ii) - L(1:jj-1,jj)' * L(1:jj-1,ii)) / L(jj,jj);
        end
    end
else
    error('Shape flag must be ''upper'' or ''lower''.');
end


