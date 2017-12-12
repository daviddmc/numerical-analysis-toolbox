function r = PolyRoot( p )
% PolyRoot  Find polynomial roots using companion matrix.
%   PolyRoot(P) computes the roots of the polynomial whose coefficients
%   are the elements of the vector P. If P has N+1 components,
%   the polynomial is P(1)*X^N + ... + P(N)*X + P(N+1).
%
%   See also PolyRootDeflation, PolyRealRootDeflation.

%   Copyright 2017 Junshen Xu   

p = p(:).';
firstNonZero = find(p, 1);
lastNonZero = find(p, 1, 'last');
if isempty(firstNonZero)
    r = [];
    return
end
p = p(firstNonZero:lastNonZero);
r = zeros(lastNonZero - firstNonZero - 1, 1);
d = p(2:end)./p(1);
while any(isinf(d))
    p = p(2:end);
    d = p(2:end)./p(1);
end

% find the roots of a polynomial as the eigenvalue of a companion matrix
n = length(p);
if n > 1
    A = diag(ones(1,n-2),-1);
    A(1,:) = -d;
    if isreal(A)
        A = DoubleHessenbergQRIter(A);
        n = size(A, 1);
        for i = 1:n-1
            if A(i+1, i) ~= 0
                a = (A(i, i) + A(i+1, i+1)) / 2;
                d = A(i, i) - A(i+1, i+1);
                bc = A(i, i+1) * A(i+1, i);
                delta = sqrt(4*bc + d^2) / 2;
                A(i,i) = a + delta;
                A(i+1, i+1) = a - delta;
            end
        end   
    else
        A = SingleHessenbergQRIter(A);
    end
    r = [r; diag(A)];
end
end

