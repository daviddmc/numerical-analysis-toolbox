function p = CauchyPoint( grad, delta, B )
%CAUCHYPOINT Summary of this function goes here
%   Detailed explanation goes here

gNorm = Norm(grad);
gBg = grad'*B*grad;
if gBg <= 0
    tau = 1;
else
    tau = min(gNorm^3 / (delta * gBg), 1);
end

p = -tau*delta/gNorm * grad;


end

