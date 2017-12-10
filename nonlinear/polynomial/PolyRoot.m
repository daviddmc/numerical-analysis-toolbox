function r = PolyRoot( p )
%POLYROOT Summary of this function goes here
%   Detailed explanation goes here

p = p(:).';
n = size(p,2);
r = [];  

inz = find(p);
if isempty(inz)
    % All elements are zero
    return
end

% Strip leading zeros and throw away.  
% Strip trailing zeros, but remember them as roots at zero.
nnz = length(inz);
p = p(inz(1):inz(nnz));
r = zeros(n-inz(nnz),1);  

% Prevent relatively small leading coefficients from introducing Inf
% by removing them.
d = p(2:end)./p(1);
while any(isinf(d))
    p = p(2:end);
    d = p(2:end)./p(1);
end

% Polynomial roots via a companion matrix
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

