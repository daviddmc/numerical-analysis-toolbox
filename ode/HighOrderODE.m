function [ output_args ] = HighOrderODE( f, a, b, y0)
%HIGHORDERODE Summary of this function goes here
%   Detailed explanation goes here

m = length(y0);


end

function dydt = F(f, x, y)
    dydt = [y(2:end), f(x, y)];
end

