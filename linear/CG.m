function x = CG( A, b ,MaxIter, x0, tolerance)
% conjugate gradient method
% solving A * x = b
% input
% A
% b
% MaxIter : max number of iteration, default: size of A
% x0 : initial value of x
% tolerance : specifies the tolerance of the method, default: 1e-6
% output
% X

    CheckSquareMatrix(A);
    if(~isvector(b))
        error('b should be a vector');
    end
    if(~iscolumn(b))
        b = b(:);
    end
    
    
    if(~exist('MaxIter', 'var') || isempty(MaxIter))
        MaxIter = size(A,1);
    end
    
    if(~exist('x0','var') || isempty(x0))
        x0 = zeros(size(A,1), 1);
    end
    if(~isvector(x0))
        error('x0 should be a vector');
    end
    if(~iscolumn(x0))
        x0 = x0(:);
    end
    CheckMultiplicationSize(A,x0, b);
    
    if(~exist('tolerance', 'var') || isempty(tolerance))
        tolerance = 1e-6;
    end

    x = x0;
    r = b - A * x;
    p = r;
    
    for k = 1 : MaxIter
        if(max(abs(r)) < tolerance)
            break;
        end
        alpha = (r'*r) ./ (p'*(A*p));
        x = x + alpha * p;
        rNew = r - alpha*A*p;
        beta = (rNew' * rNew) / (r' * r);        
        p = rNew + beta * p;
        r = rNew;
    end

end

