function [D, Q] = SymEigen( A, method )
% SymEigen    Eigenvalues and eigenvectors of symmetric matrix.
%   E = SymEigen(A) produces a column vector E containing the eigenvalues 
%   of a symmetric (Hermitian) matrix A.
% 
%   [D, Q] = SymEigen(A) produces a diagonal matrix D of eigenvalues and a 
%   unitary matrix Q so that D = Q'*A*Q, columns of Q are the corresponding
%   eigenvectors.
%
%   [D, Q] = SymEigen(A, METHOD) specifies the method to calculate
%   eigenpair. Avaliable methods are:
%       
%       QR - (default) QR iteration with double shift for real matrix and
%                      single shift for complex matrix.
%       Jacobi - Jacobi iteration. 
%       
%   See also Eigen.

%   Copyright 2017 Junshen Xu

flagQ = nargout > 1;

if ~exist('method','var')
    method = 'QR';
end

if strcmpi(method, 'QR')
    if flagQ
        [D, Q] = TridiagReduction( A ); 
        [D, Q] = TridiagQRIter(D, Q);
    else
        D = TridiagReduction(A);
        D = TridiagQRIter(D);
        D = diag(D);
    end
elseif strcmpi(method, 'Jacobi')
    if flagQ
        [D, Q] = Jacobi(A);
    else
        D = Jacobi(A);
    end
end

end