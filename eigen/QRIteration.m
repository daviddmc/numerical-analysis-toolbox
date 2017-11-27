function A = QRIteration( A , shift)
%QRITERATION Summary of this function goes here
%   Detailed explanation goes here


n = size(A, 1);
Ishift = shift * eye(n);
%{
for k = 1 : 100
    [Q, R] = QR(A - Ishift);
    A = R * Q + Ishift;
end
%}


A = HessenbergReduction(A);

for k = 1 : 20
    A = A - Ishift; 
    
    t = abs(A(1, 1))^2 + abs(A(2, 1))^2;
    if(t ~= 0)
        t = sqrt(t);
        c = A(1,1) / t;
        s = A(2,1) / t;
        signi = sign(A(1,1));
        signk = sign(A(2,1));
        if(signk == 0)
            signk = 1;
        end
        q = [c * signk s*signi*signi/signk; -s c];
        A([1, 2], :) = q * A([1, 2],:);
        A(1:3, [1, 2]) = A(1:3, [1, 2]) * q';
        %A(:,[1,2]) = A(:,[1,2]) * q'
    end
    
    for ii = 2 : n-2
       t = abs(A(ii, ii))^2 + abs(A(ii+1, ii))^2;
       if(t == 0)
           continue
       end
       t = sqrt(t);
       c = A(ii,ii) / t;
       s = A(ii+1, ii) / t;
       signi = sign(A(ii,ii));
       signk = sign(A(ii+1,ii));
       if(signk == 0)
           signk = 1;
       end
       q = [c * signk s*signi*signi/signk; -s c];
       A([ii, ii+1], ii-1:end) = q * A([ii, ii+1], ii-1:end);
       A(1:ii+2, [ii, ii+1]) = A(1:ii+2, [ii, ii+1]) * q';
       %A([ii, ii+1], :) = q * A([ii, ii+1], :);
       %A(:, [ii, ii+1]) = A(:, [ii, ii+1]) * q'
    end
    
    t = abs(A(n-1, n-1))^2 + abs(A(n, n-1))^2;
    if(t ~= 0)
        t = sqrt(t);
        c = A(n-1,n-1) / t;
        s = A(n,n-1) / t;
        signi = sign(A(n-1,n-1));
        signk = sign(A(n,n-1));
        if(signk == 0)
           signk = 1;
        end
        q = [c * signk s*signi*signi/signk; -s c];
        A([n-1, n], n-2:end) = q * A([n-1, n], n-2:end);
        %A([n-1, n], :) = q * A([n-1, n], :);
        A(:, [n-1, n]) = A(:, [n-1, n]) * q';
    end
    A = A + Ishift;
    
    for ii = 1 : n-2
        A(ii+2:end,ii) = 0;
    end
end



