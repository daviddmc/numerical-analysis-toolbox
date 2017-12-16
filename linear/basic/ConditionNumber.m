function condNum = ConditionNumber( X, type )
% condition number of square matrix
% input
% X : input matrix or vector
% type : type of norm, 1, 2, Inf, 'fro'. Default: 2.
% output
% condNum : condition number of X

CheckSquareMatrix(X, 'X');

if(~exist('type', 'var'))
    type = 2;
end

condNum = Norm(X, type) * Norm(Inverse(X), type);

