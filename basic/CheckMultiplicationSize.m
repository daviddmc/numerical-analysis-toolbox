function  CheckMultiplicationSize(A,B)
    if(size(A,2) ~= size(B, 1))
        error('multiplication size unmatch');
    end 
end

