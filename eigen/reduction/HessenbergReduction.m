function [H, Q] = HessenbergReduction(A)
% HessenbergReduction   Hessenberg form.
%   H = HessenbergRection(A) is the Hessenberg form of the matrix A. The 
%   Hessenberg form of a matrix is zero below the first subdiagonal and has
%   the same eigenvalues as A. If the matrix is symmetric (Hermitian), the 
%   form is tridiagonal, which is the same as TridiagReduction.  
%
%   [H, Q] = HessenbergReduction(A) produces a Hessenberg matrix H and a 
%   unitary matrix Q and  so that A = Q*H*Q' and Q'*Q = EYE(SIZE(Q)).
% 
%   See also BidiagReduction, TridiagReduction.

%   Copyright 2017 Junshen Xu

n = size(A, 1);
if(nargout > 1)
    Q = eye(n);
end
beta = zeros(n-2, 1);
for ii = 1 : n-2
    [v, beta(ii)] = HouseVec(A(ii+1:end, ii));
    A(ii+1:end, ii:end) = A(ii+1:end, ii:end) - beta(ii)*v*(v'*A(ii+1:end, ii:end));
    A(:,ii+1:end) = A(:,ii+1:end) - beta(ii)*(A(:,ii+1:end)*v)*v';
    A(ii+2:end, ii) = v(2:end);
end

if(nargout > 1)
    for ii = n-2:-1:1
        v = [1; A(ii+2:end, ii)];
        %Q(ii+1:end, ii+1:end) = Q(ii+1:end, ii+1:end) - beta(ii)*(Q(ii+1:end, ii+1:end)*v)*v';
        Q(ii+1:end, ii+1:end) = Q(ii+1:end, ii+1:end) - beta(ii)*v*(v'*Q(ii+1:end, ii+1:end));
    end
end

H = triu(A, -1);

function [x, beta] = HouseVec(x)
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

