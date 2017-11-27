function X = LS( A, B, method )
%LS Summary of this function goes here
%   Detailed explanation goes here

if(strcmp(method, 'normal'))
   X = CholeskySolve(A' * A, A' * B);
elseif(strcmp(method, 'QR'))
    
elseif(strcmp(method, 'SVD'))
    
else
    error('method');
end

