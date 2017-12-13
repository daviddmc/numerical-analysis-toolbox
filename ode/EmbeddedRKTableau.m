function [ A,c,b,p,FSAL,bHat,pHat] = EmbeddedRKTableau(method)
% EmbeddedRKTableau    Butcher tableau for embedded RK method.
%   [A, c, b, p,bHat, pHat] = EmbeddedRKTableau(METHOD) generates the
%   following Butcher tableau of a Embedded Runge-Kutta method:
%       
%         0|
%      c(2)|A(2,1)
%      c(3)|A(3,1)   A(3,2)
%        .
%        .
%      c(s)|A(s,1)  A(s,2)  A(s,3) ... A(s,s-1)  0
%        -----------------------------------------------
%           b(1)    b(2)    b(3)   ... b(s-1)    b(s)
%           bHat(1) bHat(2) bHat(3)... bHat(s-1) bHat(s) 
%
%   where b and bHat are corresponding to a p-order and pHat-order RK
%   method respectively. FSAL is a boolean indicating whether A,b,c satisfy
%   'first same as last' property.
%   
%   METHOD specifies the embedded RK method. The available method are:
%       
%       'DP54'  - (default) Dormand-Prince 5(4) method
%       'RKF45' - Runge-Kutta-Fehlberg 4(5) method
%       'CK45'  - Cash-Karp 4(5) method
%       'BS32'  - Bogacki-Shampine 3(2) method
%       'RKF12' - Runge-Kutta-Fehlberg 1(2) method
%       'HE12'  - Heun-Euler 1(2) method

if nargin < 1
    method = 'DP54';
end

if strcmpi(method, 'DP54')
    c = [0, 0.2, 0.3, 0.8, 8/9, 1, 1];
    A = zeros(7,7);
    A(2,1) =  0.2;
    A(3,1:2) = [3/40, 9/40];
    A(4,1:3) = [44/45, -56/15, 32/9];
    A(5,1:4) = [19372/6561, -25360/2187, 64448/6561, -212/729];
    A(6,1:5) = [9017/3168, -355/33, 46732/5247, 49/176, -5103/18656];
    A(7,1:6) = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84];
    b = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];
    bHat = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    p = 5;
    pHat = 4;
    FSAL = true;
elseif strcmpi(method, 'RKF45')
    c = [0, 0.25, 3/8, 12/13, 1, 0.5];
    A = zeros(6,6);
    A(2,1) =  0.25;
    A(3,1:2) = [3/32, 9/32];
    A(4,1:3) = [1932/2197, -7200/2197, 7296/2197];
    A(5,1:4) = [439/216, -8, 3680/513, -845/4104];
    A(6,1:5) = [-8/27, 2, -3544/2565, 1859/4104, -11/40];
    b = [25/216, 0, 1408/2565, 2197/4104, -0.2, 0];
    bHat = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55];
    p = 4;
    pHat = 5;
    FSAL = false;
elseif strcmpi(method, 'CK45')
    c = [0, 0.2, 0.3, 0.6, 1, 7/8];
    A = zeros(6,6);
    A(2,1) =  0.2;
    A(3,1:2) = [3/40, 9/40];
    A(4,1:3) = [0.3, -0.9, 1.2];
    A(5,1:4) = [-11/54, 2.5, -70/27, 35/27];
    A(6,1:5) = [1631/55296, 175/512, 575/13824, 44275/110592, 253/4096];
    bHat = [37/378, 0, 250/621, 125/594, 0, 512/1771];
    b = [2825/27648, 0, 18575/48384, 13525/55296, 277/14336, 0.25];
    p = 4;
    pHat = 5;
    FSAL = false;
elseif strcmpi(method, 'BS32')
    c = [0, 0.5, 0.75, 1];
    A = zeros(4, 4);
    A(2,1) = 0.5;
    A(3,1:2) = [0, 0.75];
    A(4,1:3) = [2/9, 1/3, 4/9];
    b = [2/9, 1/3, 4/9, 0];
    bHat = [7/24, 1/4, 1/3, 1/8];
    p = 3;
    pHat = 2;
    FSAL = true;
elseif strcmpi(method, 'RKF12')
    c = [0, 0.5, 1];
    A = zeros(3,3);
    A(2, 1) = 0.5;
    A(3,1:2) = [1/256 255/256];
    b = [1/256, 255/256, 0];
    bHat = [1/512, 255/256, 1/512];
    p = 1;
    pHat = 2;
    FSAL = true;
elseif strcmpi(method, 'HE12')
    c = [0, 1];
    A = zeros(2, 2);
    A(2, 1) = 1;
    bHat = [0.5, 0.5];
    b = [1, 0];
    p = 1;
    pHat = 2;
    FSAL = true;
else
    error('method error');
end

end

