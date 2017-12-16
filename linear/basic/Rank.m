function r = Rank( A, method )
% Rank   Matrix rank.
%   Rank(A) rovides an estimate of the number of linearly independent rows 
%   or columns of a matrix A.
%
%   Rank(A, METHOD) specifies the method. The avaliable methods are:
%
%       SVD - (default) singular value decomposition
%       QR  - QR decomposition with column pivoting
%
%   See also

%   Copyright 2017 Junshen Xu


if ~exist('method', 'var')
    method = 'svd';
end

if strcmpi(method, 'svd')
    s = SVD(A);
elseif strcmpi(method, 'qr')
    [~,R,~] = QR(A, 0);
    s = abs(diag(R));
else
    error('method should be ''svd'' or ''qr''.');
end

tol = max(size(A)) * eps(s(1));
r = sum(s > tol);

