function  CheckSquareMatrix( X , name)
    if(~ismatrix(X) || isempty(X) || size(X, 1) ~= size(X, 2))
        error([name ' should be nonempty square matrix']);
    end 
end

