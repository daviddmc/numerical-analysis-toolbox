function [S, U, V] = SVD( A, econ )
% SVD    Singular value decomposition.
%   [S, U, V] = SVD(X) produces a diagonal matrix S, of the same dimension 
%   as X and with nonnegative diagonal elements, and unitary matrices U and 
%   V so that X = U*S*V'.
%
%   S = svd(X) returns a vector containing the singular values.
%   
%   [S, U, V] = SVD(X, 0) also produces the "economy size" decomposition. 
%   If X is m-by-n with m > n, then only the first n columns of U are 
%   computed and S is n-by-n. For m < n, only the first m columns of V are 
%   computed and S is m-by-m.
%
%   See also 

%   Copyright 2017 Junshen Xu

flagU = nargout > 1;
flagV = nargout > 2;

[m, n] = size(A);
if n > m
    flagT = 1;
    A = A';
else
    flagT = 0;
end

if(~exist('econ','var'))
    econ = 1;
end

if flagU
    if flagV
        [A, U, V] = BidiagReduction(A, econ); 
        [A, U, V] = BidiagQRIter(A, U, V);
    else
        [A, U] = BidiagReduction(A, econ); 
        [A, U] = BidiagQRIter(A, U);
    end
else
    A = BidiagReduction(A, econ); 
    A = BidiagQRIter(A);
end
m = size(A, 1);

if flagV
    signA = sign(diag(A)');
    signA(signA == 0) = 1;
    V = V .* signA;
end

if ~flagU
    S = abs(diag(A));
else
    S = zeros(size(A));
    [s, idx] = sort(abs(diag(A)), 'descend');
    S(1:m+1:end) = s;
    U = U(:, idx);
    V = V(:, idx);
end

if flagT
    S = S';
    tmp = V;
    V = U;
    U = tmp;
end

end