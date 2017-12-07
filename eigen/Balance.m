function [A, D] = Balance(A, p)
% Balance   Diagonal scaling to imporve eigenvalue accuracy.
%   [B, D] = Balance(A) finds a diagonal matrix D such that B = S\A*S has, 
%   as nearly as possible, approximately equal row and column norms. T is a
%   a diagonal matrix whose elements are powers of the radix base(i.e., 2).
%   This restriction ensures there is no computational error in the
%   balancing algorithm.
%
%   [B, D] = Balance(A, p) specifies the norm to be balanced, default 2.
%
%   B = Balance(A, ...) returns the balanced matrix B.
%
%   See also

%   Reference:
%   [1]James, R., Langou, J., & Lowery, B. R. (2014). On matrix balancing 
%   and eigenvector computation. Mathematics.
%   
%   Copyright 2017 Junshen Xu

if(~exist('p','var'))
    p = 2;
end

if(~isnumeric(p) || p < 0 || isinf(p) || isnan(p))
    error('p should be number in (0, infinity).');
end

n = size(A, 1);
D = eye(n);
flag = 0;
while ~flag
    flag = 1;
    for i = 1:n
        c = Norm(A(:, i), p);
        r = Norm(A(i, :), p);
        s = c^p + r^p;
        f = 1;
        while c < r/2
            c = 2*c;
            r = r/2;
            f = 2*f;
        end
        while c>= r*2
            c = c/2;
            r = r*2;
            f = f/2;
        end
        if (c^p+r^p) < 0.95*s
            flag = 0;
            D(i,i) = f*D(i,i);
            A(:,i) = f*A(:,i);
            A(i,:) = A(i,:)/f;
        end
    end
end

end

