function E = Eigen( A )
% Eigen    Eigenvalues.
%    E = Eigen(A) produces a column vector E containing the eigenvalues of 
%    a square matrix A.

A = HessenbergReduction(A);

if isreal(A)
    A = DoubleHessenbergQRIter(A);
    n = size(A, 1);
    for i = 1:n-1
        if A(i+1, i) ~= 0
            a = (A(i, i) + A(i+1, i+1)) / 2;
            d = A(i, i) - A(i+1, i+1);
            bc = A(i, i+1) * A(i+1, i);
            r = sqrt(4*bc + d^2) / 2;
            A(i,i) = a + r;
            A(i+1, i+1) = a - r;
        end
    end   
else
    A = SingleHessenbergQRIter(A);
end

E = diag(A);


end

