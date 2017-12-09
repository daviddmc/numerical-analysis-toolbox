function [df, err] = Diff(fun, x)

% five-point midpoint formulat
% 

extrapMat = [1.000000000000000   1.000000000000000   1.000000000000000;
             1.000000000000000   0.062499987500002   0.015624995312501;
             1.000000000000000   0.003906248437500   0.000244140478516;
             1.000000000000000   0.000244140478516   0.000003814693832];
pinvMat =  [0.000529100254976  -0.041798931994960   0.499470899745022   0.541798931994962;
            -0.343684340138851  22.149550057480415  -9.649788452291601 -12.156077265049964;
            1.343154800121674 -22.107715504737246   9.149832054840520  11.614728649775055];

maxStep = 100;
r = 2.0000001;
h = max(x,0.02);
delta = maxStep*r .^(0:-1:-25);
diffCoef = [-1/(r^2 - 1), r^3/(r^2 - 1)]; 

df = zeros(size(x));
err = df;
for i = 1 : numel(x)

  fPlus = fun(x(i) + h(i) * delta);
  fMinus = fun(x(i) - h(i) * delta);
  f_del = (fPlus - fMinus)/2;
  f_del = f_del(:)';
  ne = length(delta) - 3;
  der_init = diffCoef * [f_del(1:ne); f_del(2:ne+1)];
  der_init = der_init ./ (h(i) * delta(1:ne));
  
  rhs = [der_init(1:end-3);der_init(2:end-2);der_init(3:end-1);der_init(4:end)];
  rombcoefs = pinvMat*rhs;
  der_romb = rombcoefs(1,:).';
  s = sqrt(sum((rhs - extrapMat*rombcoefs).^2));
  errors = s * 9.378218028804616;
  
  [der_romb,tags] = sort(der_romb);
  der_romb([1,2,end-1,end]) = [];
  tags([1,2,end-1,end]) = [];
  errors = errors(tags);
  [err(i),ind] = min(errors);
  df(i) = der_romb(ind);
end

end
