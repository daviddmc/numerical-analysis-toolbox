function X = TikhonovLS( A, B, lambda, L)
%TikhonovLS    Tikhonov regularized least square problem.
%   X = TikhonovLS(A, B, lambda) solves the least square problem 
%       
%       argmin_x||A*x-B||_2^2 + lambda||x||_^2
%
%   where A is a M-by-N matrix, lamabda is a non-negative real number.
%
%   X = TikhonovLS(A, B, lambda, L) solves the least square problem 
%       
%       argmin_x||A*x-B||_2^2 + lambda||L*x||_^2
%
%   See also

%   Copyright 2017 Junshen Xu

if ~exist('L','var')
    L = eye(size(A, 2));
end

if ~isreal(lambda) || lambda < 0
    error('lambda should be non-negative real number');
end

A = [A; sqrt(lambda) * L];
B = [B; zeros(size(L, 1), size(B,2))];

X = LS(A, B);

end

