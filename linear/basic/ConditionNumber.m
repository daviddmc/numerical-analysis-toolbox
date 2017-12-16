function condNum = ConditionNumber( X, type )
% ConditionNumber   Condition number with respect to inversion.
%   cond(X) returns the 2-norm condition number (the ratio of the largest 
%   singular value of X to the smallest). 
%   
%   cond(X,P) returns the condition number of X in P-norm (matrix norm).
%
%   See also

%   Copyright 2017 Junshen Xu

CheckSquareMatrix(X, 'X');

if(~exist('type', 'var'))
    type = 2;
end

if type == 2
    s = SVD(X);
    condNum = s(1) / s(end);
else
    condNum = Norm(X, type) * Norm(Inverse(X), type);
end

