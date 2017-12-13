function [ A,c,b,p] = RKTableau(method)
% RKTableau    Butcher tableau for RK method.
%   [A, c, b, p] = RKTableau(METHOD) generates the following Butcher 
%   tableau of a Embedded Runge-Kutta method:
%       
%         0|
%      c(2)|A(2,1)
%      c(3)|A(3,1)   A(3,2)
%        .
%        .
%      c(s)|A(s,1)  A(s,2)  A(s,3) ... A(s,s-1)  0
%        -----------------------------------------------
%           b(1)    b(2)    b(3)   ... b(s-1)    b(s)
%
%   where b is corresponding to a p-order RK method.
%   
%   METHOD specifies the embedded RK method. The available method are:
%       
%       'Euler'    - Euler method (1-order)
%       'midpoint' - midpoint method (2-order)
%       'Heun'     - Heun method (2-order)
%       'Ralston'  - Ralston method (2-order)
%       'Kutta3'   - 3-order Kutta method
%       'classic4' - (default) classic 4-order method

if nargin < 1
    method = 'classic4';
end

if strcmpi(method, 'Euler')
    c = 0;
    A = 0;
    b = 1;
    p = 1;
elseif strcmpi(method, 'midpoint')
    c = [0, 0.5];
    A = [0 0; 0.5 0];
    b = [0, 1];
    p = 2;
elseif strcmpi(method, 'Heun')
    c = [0, 1];
    A = [0, 0; 1, 0];
    b = [0.5, 0.5];
    p = 2;
elseif strcmpi(method, 'Ralston')
    c = [0, 2/3];
    A = [0, 0; 2/3, 0];
    b = [0.25, 0.75];
    p = 2;
elseif strcmpi(method, 'Kutta3')
    c = [0, 0.5, 1];
    A = zeros(3);
    A(2,1) = 0.5;
    A(3,1:2) = [-1, 2];
    b = [1/6, 2/3, 1/6];
    p = 3;
elseif strcmpi(method, 'classic4')
    c = [0, 0.5, 0.5, 1];
    A = sparse(diag([0.5, 0.5, 1], -1));
    b = [1/6, 1/3, 1/3, 1/6];
    p = 4;
else
    error('method error');
end

end

