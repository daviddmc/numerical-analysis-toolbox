function normX = Norm( X, type )
% Norm   Matrix or vector norm.
%   norm(X,2) returns the 2-norm of X.
%
%   norm(X) is the same as norm(X,2).
%
%   norm(X,1) returns the 1-norm of X.
%   
%   norm(X,Inf) returns the infinity norm of X.
%   
%   norm(X,'fro') returns the Frobenius norm of X.
%   
%   In addition, for vectors...
%   
%   norm(V,P) returns the p-norm of V defined as SUM(ABS(V).^P)^(1/P).
%   
%   norm(V,Inf) returns the largest element of ABS(V).
%   
%   norm(V,-Inf) returns the smallest element of ABS(V).
%
%   See also

%   Copyright 2017 Junshen Xu

if(~ismatrix(X))
    error('input X should be vector or matrix');
end

if(isempty(X))
    normX = 0;
    return
end

if ~exist('type', 'var')
    type = 2;
end

if isvector(X)
    if(type == 1)
        normX = sum(abs(X));
    elseif(type == 2 || strcmpi(type, 'fro'))
        normX = sqrt(sum(abs(X).^2));
    elseif(type == inf)
        normX = max(abs(X));
    elseif(type == -inf)
        normX = min(abs(X));
    elseif(isnumeric(type))
        normX = (sum(abs(X))^type)^(1/type);
    else
        error('type of vector norm should be real number, inf, -inf or ''fro''.');
    end
else
    if(type == 1)
        normX = max(sum(abs(X), 1));
    elseif(type == 2)
        s = SVD(X);
        normX = s(1);
    elseif(type == inf)
        normX = max(sum(abs(X), 2));
    elseif(type == 'fro')
        normX = sqrt(sum(sum(abs(X).^2)));
    else
        error('type of matrix norm should be 1, 2, inf or ''fro''.');
    end
end




