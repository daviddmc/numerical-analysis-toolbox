function [tout, yout] = RK(OdeFcn,tSpan,y)

% ------- INPUT PARAMETERS
% Time properties
tSpan  = tSpan(:);                     % Column vector
t      = tSpan(1);
tEnd = tSpan(end);
tDir = sign(tEnd-t);

% Number of equations, y is a column vector
Ny     = length(y);

% General options
AbsTol      = 1e-6;
RelTol      = 1e-3;
%h           = 1e-2;         % h may be positive or negative
hmax        = abs(tSpan(end) - tSpan(1));        % hmax is positive
MaxNbrStep  = 1e5;
FacL        = 1/0.2;
FacR        = 1/10;
Safe        = 0.9;%0.9;                % WORK(2)
NormControl = false;
Beta        = 0;%0.04;

nBuffer = 100;
nout    = 0;
tout    = zeros(nBuffer,1);
yout    = zeros(nBuffer,Ny);
 

% ---------------------------
% Coefficients values for dop54
% ---------------------------
[a,c,b,p] = RKTableau('classic4');
s = length(c);

% Initialiation of internal parameters
FacOld = 1e-4;
Expo1  = 1/(p + 1) - Beta*0.75;             
K      = zeros(Ny,s);  % For DOP54 a is a matrix(7,7)
% --------------
% Integration step
f0 = feval(OdeFcn,t,y);
K(:,1) = f0;
hmax = min(hmax,abs(tSpan(end)-tSpan(1)));      % hmax positive
% init step
h = hInitFcn(OdeFcn,t,y,tDir,hmax,K(:,1),RelTol,AbsTol, p);
h = tDir * min(abs(h),hmax);
Reject = false;

% ------------
% --- BASIC INTEGRATION STEP  
% ------------

Done = false;
StepAccNbr = 0;
for iter = 1:MaxNbrStep
  if Done
      break
  end
  
  if (0.1*abs(h) <= abs(t)*eps)
    warning('Too small step size');    
    Done = true;
  end
  
  if  (tDir * (t + 1.001*h - tEnd) >= 0 )     
    h = tEnd - t;
  elseif tDir * (t + 1.8*h - tEnd) > 0
    h = (tEnd - t)*0.5;
  end
  
  ch = h*c;  
  ah = h*a';     % Needed for matrix calculation
  bh = h*b;
  
  for j = 2:s  
    y1 = y + K * ah(:, j);
    f0 = feval(OdeFcn,t+ch(j), y1);
    K(:,j) = f0; 
  end
  yNew = y + K*bh';
  
  ahHalf = ah / 2;
  chHalf = ch / 2;
  bhHalf = bh / 2;
  for j = 2:s  
    yHalf1 = y + K * ahHalf(:, j);
    f0 = feval(OdeFcn,t+chHalf(j), yHalf1);
    K(:,j) = f0; 
  end
  yHalf1 = y + K * bhHalf';
  
  tphh = t + h /2;
  K(:, 1) = feval(OdeFcn,tphh, yHalf1);

  for j = 2:s
      yHalf2 = yHalf1 + K * ahHalf(:, j);
      f0 = feval(OdeFcn,tphh+chHalf(j), yHalf2);
      K(:,j) = f0;
  end

  yHalf2 = yHalf1 + K * bhHalf';
  K4 = (yHalf2 - yNew) / (2^p - 1);
  yNew = yHalf2 + K4;
  
  tph  = t + h;
  % --- ERROR ESTIMATION   (450)
  if NormControl
    %  norm(e) <= max(RelTol*norm(y),AbsTol)
    Sk = max(AbsTol) + max(RelTol)* max(norm(y),norm(yNew));
    Err = norm(K4)/Sk;
  else  
    Sk   = AbsTol + RelTol .* max( abs(y),abs(yNew));         
    Err  = sqrt( sum((K4./ Sk).^2)/Ny );
  end
  % --- COMPUTATION OF HNEW -----> 662 Hairer
  Fac11 = Err^Expo1;
  % --- LUND-STABILIZATION
  Fac   = Fac11/FacOld^Beta;
  % --- WE REQUIRE  FAC1 <= HNEW/H <= FAC2
  Fac   = max(FacR,min(FacL,Fac/Safe));
  hNew  = h/Fac;   
  if(Err < 1.D0)
    % --- STEP IS ACCEPTED          (470 Hairer)
    StepAccNbr = StepAccNbr + 1;
    FacOld = max(Err,1e-4);
    % ------- STIFFNESS DETECTION                     675
    
    
% ------- FINAL PREPARATION FOR DENSE OUTPUT     495         
    nout = nout + 1;
    if nout > length(tout)
        tout = [tout;zeros(nBuffer,1)];
        yout = [yout;zeros(nBuffer,Ny)];  
    end
    tout(nout)   = t;
    yout(nout,:) = y';
    
    t = tph;
    y = yNew;  
    
    if t == tEnd  
      Done = true;
    end           
    
    if abs(hNew) > hmax
      hNew = tDir*hmax;
    end
    if Reject
      hNew = tDir*min(abs(hNew),abs(h));
      Reject = false;
    end	     
    
  else % --- STEP IS REJECTED      depuis 457    (769 Hairer)  
    hNew = h/min(FacL,Fac11/Safe);
    Reject = true;    
  end
  h = hNew;
  K(:,1) = feval(OdeFcn,t,y);
end  % main loop

nout = nout + 1;
tout(nout) = t;
yout(nout,:) = y';
tout = tout(1:nout);
yout = yout(1:nout,:);



function h = hInitFcn(OdeFcn,t,y,PosNeg,hmax,f0,RelTol,AbsTol,order)                          

Sk  = AbsTol + RelTol.*abs(y);
Dnf = sum( (f0./Sk).^2 );
Dny = sum( (y./Sk).^2 );
if (Dnf < 1e-10 || Dny < 1e-10)
  h = 1e-6;
else
  h = sqrt(Dny/Dnf) * 0.01; 
end
h = min(h,hmax) * PosNeg;

% ---- PERFORM AN EXPLICIT EULER STEP
y1 = y + h*f0;
f1  = feval(OdeFcn,t+h,y1); 

% ---- ESTIMATE THE SECOND DERIVATIVE OF THE SOLUTION
Sk   = AbsTol + RelTol .* abs(y);
Der2 = sum ( ((f1-f0)./Sk).^2 );   
Der2 = sqrt(Der2)/h;
% ---- STEP SIZE IS COMPUTED SUCH THAT
% ----  H**IORD * MAX ( NORM (F0), NORM (DER2)) = 0.01
Der12 = max(abs(Der2),sqrt(Dnf));
if Der12 <= 1e-15 
  h1 = max(1e-6,abs(h)*1e-3);
else
  h1 = (0.01/Der12)^(1/order);
end
h = min([100*abs(h),h1,hmax])*PosNeg;
return

 
