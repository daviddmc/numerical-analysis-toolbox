function X = TikhonovLS( A, B, lambda, L)


if ~exist('L','var')
    L = eye(size(A, 2));
end

A = [A; sqrt(lambda) * L];
B = [B; zeros(size(L, 1), size(B,2))];

X = LS(A, B);


end

