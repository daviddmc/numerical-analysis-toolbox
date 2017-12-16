function normX = Norm( X, type )
% Matrix or vector norm
% input
% X : input matrix or vector
% type : type of norm, 1, 2, Inf, 'fro'. Default: 2.
% output
% normX : norm of X

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

if(type == 1)
    if(isvector(X))
        normX = sum(abs(X));
    else
        normX = max(sum(abs(X), 1));
    end
elseif(type == 2)
    if(isvector(X))
        normX = sqrt(sum(abs(X).^2));
    else
        s = SVD(X);
        normX = s(1);
    end
elseif(type == inf)
    if(isvector(X))
        normX = max(abs(X));
    else
        normX = max(sum(abs(X), 2));
    end
elseif(type == 'fro')
    normX = sqrt(sum(sum(abs(X).^2)));
elseif(isnumeric(type))
    if isvector(X)
        normX = (sum(abs(X))^type)^(1/type);
    else
        error('The only matrix norms available are 1, 2, inf, and ''fro''.');
    end
else
    error('type of norm should be 1, 2, inf or ''fro''.');
end

