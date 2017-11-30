function x = FFT( x )
%FFT Summary of this function goes here
%   Detailed explanation goes here

x = Base3FFT(x);

end

function x = Base2FFT(x)


flag = 0;
if(iscolumn(x))
    flag = 1;
    x = x.';
end

N = length(x);
if mod(N, 2) > 0 
    error(' ');
end

N_min = min(N, 128);
n = 0 : N_min-1;
k = n';
M = exp(-2j * pi * bsxfun(@times, n, k) / N_min);
x = reshape(x, [N/N_min, N_min]) * M;
%{
while(size(x,1) < N)
    x_even = x(:, 1:size(x, 2) /2);
    x_odd = x(:, size(x,2)/2+1:end);
    factor = exp(-1j * pi * (0:size(x,1)-1)' / size(x,1));
    x = [x_even + factor .* x_odd; x_even - factor .* x_odd];    
end
%}
while(size(x,2) < N)
    x_even = x(1:size(x, 1) /2, :);
    x_odd = x(size(x,1)/2+1:end, :);
    factor = exp(-1j * pi * (0:size(x,2)-1) / size(x,2));
    x = [x_even + factor .* x_odd, x_even - factor .* x_odd];    
end

if(flag)
    x = x.';
end

end

function x = Base3FFT(x)

if(iscolumn(x))
    x = x';
end

N = length(x);
if mod(N, 3) > 0 
    error(' ');
end

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

if(flag)
    x = x.';
end

end