function  CheckMultiplicationSize(A, X, B)
    
    if(~isempty(X))
        if(size(X,1) ~= size(A,2))
            error('size unmatch');
        end
        
        if(~isempty(B))
            if(size(X,2) ~= size(B,2))
                error('size unmatch');
            end
        end
    end
    
    if(~isempty(B))
        if(size(B,1) ~= size(A,1))
            error('size unmatch');
        end
    end
end

