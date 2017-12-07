function [B, U, V] = BidiagReduction(A, econ)
% BidiagReduction   Tridiagonal form.
%   B = BidiagReduction(A) is the bidiagonal form of the M-by-N matrix A.
%   When M >= N, T is a upper triangle and is zero above the first
%   superdiagonal. Otherwise T is a lower triangle whose elements below the
%   first subdiagonal are zero.
%
%   [B, U, V] = BidiagReduction(A) produces a bidiagonal matrix T and a 
%   unitary matrix U and V so that U'*A*V = B.
%
%   [B, U, V] = BidiagReduction(A, 0) produces the "economy size" 
%   decomposition.
% 
%   See also TridiagReduction, HessenbergReduction.

%   Copyright 2017 Junshen Xu

[m, n] = size(A);

flagU = nargout > 1;
flagV = nargout > 2;

if(nargin > 1 && econ == 0)
    flagEcon = 1;
else
    flagEcon = 0;
end

if m >= n
    beta1 = zeros(n, 1);
    beta2 = zeros(n-2, 1);
    for j = 1:n
        [v, beta1(j)] = HouseVec(A(j:end, j));
        A(j:end, j:end) = A(j:end, j:end) - beta1(j)*v*(v'*A(j:end, j:end));
        if flagU
            A(j+1:end, j) = v(2:end);
        end
        if j <= n-2
            [v, beta2(j)] = HouseVec(A(j,j+1:end)');
            A(j:end, j+1:end) = A(j:end, j+1:end) - beta2(j) *(A(j:end, j+1:end)*v)*v';
            if flagV
                A(j,j+2:end) = v(2:end)';
            end
        end
    end
else
    beta1 = zeros(m-2, 1);
    beta2 = zeros(m, 1);
    for j = 1:m
        [v, beta2(j)] = HouseVec(A(j,j:end)');
        A(j:end, j:end) = A(j:end, j:end) - beta2(j)*(A(j:end, j:end)*v)*v';
        if flagV
            A(j, j+1:end) = v(2:end)';
        end
        if j <= m-2
            [v, beta1(j)] = HouseVec(A(j+1:end, j));
            A(j+1:end, j:end) = A(j+1:end, j:end) - beta1(j) *v*(v'*A(j+1:end,j:end));
            if flagU
                A(j+2:end, j) = v(2:end);
            end
        end
    end
end

if flagU
    if flagEcon
        U = eye(m, min(m,n));
    else
        U = eye(m);
    end
    k = m < n;
    for j = length(beta1) : -1 : 1
        v = [1; A(j+1+k:end, j)];
        U(j+k:end,j+k:end) = U(j+k:end,j+k:end) - beta1(j)*v*(v'*U(j+k:end,j+k:end));
    end
end

if flagV
    if flagEcon
        V = eye(n, min(m,n));
    else
        V = eye(n);
    end
    k = m >= n;
    for j = length(beta2):-1:1
        v = [1; A(j,j+1+k:end)'];
        V(j+k:end,j+k:end) = V(j+k:end,j+k:end) - beta2(j)*v*(v'*V(j+k:end,j+k:end));
    end
end

if flagEcon
    B = zeros(min(m,n));
    B(1:size(B,1)+1:end) = diag(A);
    if m >= n
        B(size(B,1)+1:size(B,1)+1:end) = A(m+1:m+1:end);
    else
        B(2:size(B,1)+1:end) = A(2:m+1:m^2);
    end
else
    B = zeros(size(A));
    if m >= n
        B(1:m+1:end) = diag(A);
        B(m+1:m+1:end) = A(m+1:m+1:end);
    else
        B(1:m+1:m^2) = diag(A);
        B(2:m+1:m^2) = A(2:m+1:m^2);
    end
end


end

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
end

