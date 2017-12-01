function x = FFT( x )
%FFT Summary of this function goes here
%   Detailed explanation goes here

N = length(x);

if N <= 128
    n = 0 : N-1;
    k = n';
    M = exp(-2j * pi * bsxfun(@times, n, k) / N);
    x = reshape(x, [1, N]) * M;
else
    [N1, N2, flag] = FactorHelper(N);
        
    x = reshape(x, [N1, N2]);
    
    if(N2 > 1)
        for ii = 1 : N1
            x(ii, :) = FFT(x(ii, :));
        end
        fac = exp(-2j * pi / N * bsxfun(@times, (0: N1-1)',0 : N2-1));
        x = fac .* x;
    end
    
    switch(flag)
        case 0
            for jj = 1 : N2
                x(:, jj) = FFT(x(:, jj));
            end
        case 1
            for jj = 1 : N2
                x(:, jj) = Base2FFT(x(:, jj));
            end
        case 2
            for jj = 1 : N2
                x(:, jj) = Base3FFT(x(:, jj));
            end
    end

    x = x.';
    x = x(:);
end

end

function x = Base2FFT(x)

N = length(x);

N_min = min(N, 128);
n = 0 : N_min-1;
k = n';
M = exp(-2j * pi * bsxfun(@times, n, k) / N_min);
x = reshape(x, [N/N_min, N_min]) * M;

while(size(x,2) < N)
    x_even = x(1:size(x, 1) /2, :);
    x_odd = x(size(x,1)/2+1:end, :);
    factor = exp(-1j * pi * (0:size(x,2)-1) / size(x,2));
    x = [x_even + factor .* x_odd, x_even - factor .* x_odd];    
end

end

function x = Base3FFT(x)

N = length(x);

N_min = min(N, 81);
n = 0 : N_min-1;
k = n';
M = exp(-2j * pi * bsxfun(@times, n, k) / N_min);
x = reshape(x, [N/N_min, N_min]) * M;

while(size(x,2) < N)
    x1 = x(1:size(x, 1) /3,:);
    x2 = x(size(x, 1)/3+1 : 2*size(x,1)/3, :);
    x3 = x(2*size(x, 1)/3+1 : end, :);
    factor1 = exp(-2j / 3 * pi * (0:size(x,2)-1) / size(x,2));
    factor2 = factor1.^2;
    x = [x1 +                 factor1 .* x2 +                 factor2 .* x3, ...
         x1 + exp(-2j/3*pi) * factor1 .* x2 + exp(-4j/3*pi) * factor2 .* x3,...
         x1 + exp(-4j/3*pi) * factor1 .* x2 + exp(-2j/3*pi) * factor2 .* x3];    
end

end

function [N1, N2, flag] = FactorHelper(N)
flag = 0;
f = factor(N);
e2 = find(f == 2);
if(~isempty(e2))
    e2 = e2(end) - e2(1) + 1;
    N1 = 2^e2;
    if (N1 >= 8)
        flag = 1;
    end   
end

if(~flag)
    e3 = find(f == 3);
    if(~isempty(e3))
        e3 = e3(end) - e3(1) + 1;
        N1 = 3^e3;
        if(N1 >= 9)
            flag = 2;
        end
    end
end

if(~flag)
    N1 = 1;
    for ii = 1:length(f)
        N1 = N1 * f(ii);
        if(N1^2 > N)
            break;
        end
    end
end
N2 = N / N1;
end