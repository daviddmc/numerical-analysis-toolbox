function X = TridiagSolve(a, b, c, D)
% TridiagSolve     Solve a tridiagonal linear system with Thomas method.
%   X = TridiagSolve(A, B, C, D) solves a tridiagonal linear system with
%   Thomas method. A, B, C are three vectors represent the sub-, main and 
%   super- diagonal of the coefficient matrix. Let the lenght of B be N. If
%   if the length of A >= N, then the following system will be solved.
%
%   [B1 C1           ]
%   [A2 B2 C2        ]
%   [   A3 B3 \\     ]  * X = D
%   [      \\ \\ CN-1]
%   [        AN    BN]
%
%   If A is of length N-1, then the following system will be solved.
%
%   [B1 C1           ]
%   [A1 B2 C2        ]
%   [   A2 B3 \\     ]  * X = D
%   [      \\ \\ CN-1]
%   [        AN-1  BN]
%
%   Note that this algorithm is implemented without pivoting.
%
%   See also GaussElimination, SPDSolver, TriangleSolve.

%   Copyright 2017 Junshen Xu

n = length(b);
if(length(a) < n-1 || length(c) < n-1)
    error('length(A) and length(C) should >= length(B) - 1');
end

if(size(D, 1) ~= n)
    error('rhs size mismatching');
end

u = zeros(n , 1);
l = zeros(n , 1);
X = D;
u(1) = b(1);
if(length(a) >= n)
    for ii = 2:n
        l(ii) = a(ii) / u(ii - 1);
        u(ii) = b(ii) - l(ii) * c(ii - 1);
    end
else
    for ii = 2:n
        l(ii) = a(ii - 1) / u(ii - 1);
        u(ii) = b(ii) - l(ii) * c(ii - 1);
    end
end

for ii = 2 : n
    X(ii, :) = X(ii,:) - l(ii) * X(ii-1,:);
end

X(n,:) = X(n, :) ./ u(n);
for ii = n-1 : -1 : 1
    X(ii,:) = (X(ii,:) - c(ii) * X(ii+1,:))/u(ii);
end






