function [D, Q] = SymEigen( A )
% SymEigen    Eigenvalues and eigenvectors of symmetric matrix.
%   E = SymEigen(A) produces a column vector E containing the eigenvalues 
%   of a symmetric (Hermitian) matrix A.
% 
%   [D, Q] = SymEigen(A) produces a diagonal matrix D of eigenvalues and a 
%   unitary matrix Q so that D = Q'*A*Q, columns of Q are the corresponding
%   eigenvectors.

flagQ = nargout > 1;

if flagQ
    [D, Q] = TridiagReduction( A ); 
    [D, Q] = TridiagQRIter(D, Q);
else
    D = TridiagReduction(A);
    D = TridiagQRIter(D);
    D = diag(D);
end

end