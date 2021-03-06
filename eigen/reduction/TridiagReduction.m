function [T, Q] = TridiagReduction( A )
% TridiagReduction   Symmetric tridiagonal form.
%   T = TridiagReduction(A) is the tridiagonal form of the symmetric
%   (Hermitian) matrix A. This function produces the same results when the
%   input matrix is symmetric (Hermitian), but faster.
%
%   [H, Q] = TridiagReduction(A) produces a tridiagonal matrix T and a 
%   unitary matrix Q and  so that A = Q*T*Q' and Q'*Q = EYE(SIZE(Q)).
% 
%   See also BidiagReduction, HessenbergReduction.

%   Copyright 2017 Junshen Xu

n = size(A,1);
beta = zeros(n-2,1);
for k = 1:n-2
    [v, beta(k), mu] = HouseVec(A(k+1:end,k));
    p = beta(k)*A(k+1:n,k+1:n)*v;
    w = p - (beta(k)*(p'*v)/2)*v;
    A(k+1, k) = mu;
    %A(k,k+1) = A(k+1,k)';
    A(k+1:n,k+1:n) = A(k+1:n,k+1:n) - v*w'-w*v';
    A(k+2:n,k) = v(2:end); 
end

if nargout > 1
    Q = eye(n);
    for k = n-2:-1:1
        v = [1; A(k+2:n,k)];
        Q(k+1:end, k+1:end) = Q(k+1:end, k+1:end) - beta(k)*v*(v'*Q(k+1:end, k+1:end));
    end
end

T = diag(diag(A));
T(2:n+1:end) = A(2:n+1:end);
T(n+1:n+1:end) = conj(T(2:n+1:end));


function [x, beta, mu] = HouseVec(x)
sigma = Norm(x);
if(sigma == 0)
    beta = 0;
    return;
else
    beta = 1 / (sigma*(sigma + abs(x(1))));
end
    
if x(1) == 0
    mu = sigma;
else
    mu = -x(1)/abs(x(1)) * sigma;
end

x(1) = x(1) - mu;
beta = beta * abs(x(1))^2;
x = x / x(1);