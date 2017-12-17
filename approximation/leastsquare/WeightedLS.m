function X = WeightedLS( A, B, D )
%WeightedLS    Weighted least square problem.
%   X = WeightedLS(A, B, D) solves the least square problem 
%       
%       argmin_x D||A*x-B||_2^2
%
%   where A is a M-by-N matrix, D is a vector of length M.
%
%   See also

%   Copyright 2017 Junshen Xu

if ~isreal(D) || any(D < 0)
    error('D should be a real non-negative vector')
end

D = sqrt(D(:));
A = D.*A;
B = D.*B;

X = LS(A, B);

end

