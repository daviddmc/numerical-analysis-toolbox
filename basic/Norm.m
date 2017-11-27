function normX = Norm( X, type )
% Matrix or vector norm
% input
% X : input matrix or vector
% type : type of norm, 1, 2, Inf, 'fro'. Default: 2.
% output
% normX : norm of X

if(~ismatrix(X) || isempty(X))
    error('input X should be nonempty vector or matrix');
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
        normX = norm(X, 2);
    end
elseif(type == inf)
    if(isvector(X))
        normX = max(abs(X));
    else
        normX = max(sum(abs(X), 2));
    end
elseif(type == 'fro')
    normX = sqrt(sum(sum(abs(X).^2)));
else
    error('type of norm should be 1, 2, inf or ''fro''.');
end

