function [tout, yout] = adams(odeFcn,tspan,y0)

neq = length(y0);
t0 = tspan(1); 
tfinal = tspan(end); 
tdir = sign(tfinal - t0);
f0 = feval(odeFcn,t0,y0); 
htspan = abs(tspan(2) - tspan(1));
hmax = abs(0.1*(tfinal-t0));

atol = 1e-6;
rtol = 1e-3;
threshold = atol / rtol;

% Handle the output

refine = 1;
t = t0;
y = y0;
yp = f0;  

% Allocate memory if we're generating output.

chunk = min(max(100,50*refine), refine+floor((2^13)/neq));
tout = zeros(1,chunk);
yout = zeros(neq,chunk);

nout = 1;
tout(nout) = t;
yout(:,nout) = y;  

% Initialize method parameters.
maxk = 12;
two = 2 .^ (1:13)';
gstar = [ 0.5000;  0.0833;  0.0417;  0.0264;  ...
          0.0188;  0.0143;  0.0114;  0.00936; ...
          0.00789;  0.00679; 0.00592; 0.00524; 0.00468];

hmin = 16*eps(t);

% Compute an initial step size h using y'(t).
absh = min(hmax, htspan);
rh = norm(yp ./ max(abs(y),threshold),inf) / (0.25 * sqrt(rtol));
if absh * rh > 1
    absh = 1 / rh;
end
absh = max(absh, hmin);

% Initialize.
k = 1;
K = 1;
phi = zeros(neq,14);
phi(:,1) = yp;
psi = zeros(12,1);
alpha = zeros(12,1);
beta = zeros(12,1);
sig = zeros(13,1);
sig(1) = 1;
w = zeros(12,1);
v = zeros(12,1);
g = zeros(13,1);
g(1) = 1;
g(2) = 0.5;

hlast = 0;
klast = 0;
phase1 = true;

% THE MAIN LOOP

done = false;
while ~done
  
  % By default, hmin is a small number such that t+hmin is only slightly
  % different than t.  It might be 0 if t is 0.
  hmin = 16*eps(t);
  absh = min(hmax, max(hmin, absh));    % couldn't limit absh until new hmin
  h = tdir * absh;
  
  % Stretch the step if within 10% of tfinal-t.
  if 1.1*absh >= abs(tfinal - t)
    h = tfinal - t;
    absh = abs(h);
    done = true;
  end
   
  % LOOP FOR ADVANCING ONE STEP.
  failed = 0;
  invwt = 1 ./ max(abs(y),threshold);
  while true
    % Compute coefficients of formulas for this step.  Avoid computing
    % those quantities not changed when step size is not changed.

    % ns is the number of steps taken with h, including the 
    % current one.  When k < ns, no coefficients change
    if h ~= hlast  
      ns = 0;
    end
    if ns <= klast 
      ns = ns + 1;
    end
    if k >= ns
      beta(ns) = 1;
      alpha(ns) = 1 / ns;
      temp1 = h * ns;
      sig(ns+1) = 1;
      for i = ns+1:k
        temp2 = psi(i-1);
        psi(i-1) = temp1;
        temp1 = temp2 + h;
        beta(i) = beta(i-1) * psi(i-1) / temp2;
        alpha(i) = h / temp1;
        sig(i+1) = i * alpha(i) * sig(i);
      end
      psi(k) = temp1;

      % Compute coefficients g.
      if ns == 1                        % Initialize v and set w
        v = 1 ./ (K .* (K + 1));
        w = v;
      else
        % If order was raised, update diagonal part of v.
        if k > klast
          v(k) = 1 / (k * (k+1));
          for j = 1:ns-2
            v(k-j) = v(k-j) - alpha(j+1) * v(k-j+1);
          end
        end
        % Update v and set w.
        for iq = 1:k+1-ns
          v(iq) = v(iq) - alpha(ns) * v(iq+1);
          w(iq) = v(iq);
        end
        g(ns+1) = w(1);
      end

      % Compute g in the work vector w.
      for i = ns+2:k+1
        for iq = 1:k+2-i
          w(iq) = w(iq) - alpha(i-1) * w(iq+1);
        end
        g(i) = w(1);
      end
    end   

    % Change phi to phi star.
    i = ns+1:k;
    phi(:,i) = phi(:,i) * diag(beta(i));

    % Predict solution and differences.
    phi(:,k+2) = phi(:,k+1);
    phi(:,k+1) = zeros(neq,1);
    p = zeros(neq,1);
    for i = k:-1:1
        p = p + g(i) * phi(:,i);
        phi(:,i) = phi(:,i) + phi(:,i+1);
    end
    
    p = y + h * p;
    tlast = t;
    t = tlast + h;
    if done
      t = tfinal;   % Hit end point exactly.
    end
    yp = feval(odeFcn,t,p);

    % Estimate errors at orders k, k-1, k-2.
    phikp1 = yp - phi(:,1);
    temp3 = norm(phikp1 .* invwt,inf);
    err = absh * (g(k) - g(k+1)) * temp3;
    erk = absh * sig(k+1) * gstar(k) * temp3;
    if k >= 2
        erkm1 = absh*sig(k)*gstar(k-1)*norm((phi(:,k)+phikp1) .* invwt,inf);
    else
        erkm1 = 0.0;
    end
    if k >= 3
        erkm2 = absh*sig(k-1)*gstar(k-2)*norm((phi(:,k-1)+phikp1) .* invwt,inf);
    else
        erkm2 = 0.0;
    end
    
    % Test if order should be lowered
    if (k == 2) && (erkm1 <= 0.5*erk)
      knew = k - 1;
    elseif (k > 2) && (max(erkm1,erkm2) <= erk)
      knew = k - 1;
    else
        knew = k;
    end
    
    % Test if step successful
    if err > rtol                       % Failed step           
      if absh <= hmin
        warning('too small step size');  
        return;
      end
      
      % Restore t, phi, and psi.
      phase1 = false;
      t = tlast;
      for i = K
        phi(:,i) = (phi(:,i) - phi(:,i+1)) / beta(i);
      end
      for i = 2:k
        psi(i-1) = psi(i) - h;
      end

      failed = failed + 1;
      reduce = 0.5;
      if failed == 3
        knew = 1;
      elseif failed > 3
        reduce = min(0.5, sqrt(0.5*rtol/erk));
      end
      absh = max(reduce * absh, hmin);
      h = tdir * absh;
      k = knew;
      K = 1:k;
      done = false;
    else                                % Successful step
      break;
    end
  end
                   
  klast = k;
  hlast = h;

  % Correct and evaluate.
  y = p + h * g(k+1) * phikp1;
  yp = feval(odeFcn,t,y);              
  
  % Update differences for next step.
  phi(:,k+1) = yp - phi(:,1);
  phi(:,k+2) = phi(:,k+1) - phi(:,k+2);
  for i = K
    phi(:,i) = phi(:,i) + phi(:,k+1);
  end

  if (knew == k-1) || (k == maxk)
    phase1 = false;
  end

  % Select a new order.
  kold = k;
  if phase1                             % Always raise the order in phase1
    k = k + 1;
  elseif knew == k-1                    % Already decided to lower the order
    k = k - 1;
    erk = erkm1;
  elseif k+1 <= ns                      % Estimate error at higher order
    erkp1 = absh * gstar(k+1) * norm(phi(:,k+2) .* invwt,inf);
    if k == 1
      if erkp1 < 0.5*erk
        k = k + 1;
        erk = erkp1;
      end
    else
      if erkm1 <= min(erk,erkp1)
        k = k - 1;
        erk = erkm1;
      elseif (k < maxk) && (erkp1 < erk)
        k = k + 1;
        erk = erkp1;
      end
    end
  end
  if k ~= kold
    K = 1:k;
  end

  nout = nout + 1;
  if nout > length(tout)
    tout = [tout, zeros(1,chunk)];  % requires chunk >= refine
    yout = [yout, zeros(neq,chunk)];
  end       
  tout(nout) = t;
  yout(:,nout) = y;
  
  if done
    break
  end
  
  % Select a new step size.
  if phase1
    absh = 2 * absh;
  elseif 0.5*rtol >= erk*two(k+1)
    absh = 2 * absh;      
  elseif 0.5*rtol < erk
    reduce = (0.5 * rtol / erk)^(1 / (k+1));
    absh = absh * max(0.5, min(0.9, reduce));
  end
  
end

tout = tout(1:nout);
yout = yout(:, 1:nout);