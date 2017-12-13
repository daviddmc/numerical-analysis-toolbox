function [dFdy,fac] = Jacobian(F,y,Fy,thresh,fac,vectorized)


% Initialize.
Fy = Fy(:);
br = eps ^ (0.875);
bl = eps ^ (0.75);
bu = eps ^ (0.25);
facmin = eps ^ (0.78);
facmax = 0.1;
ny = length(y);
nF = length(Fy);
if isempty(fac)
  fac = sqrt(eps) + zeros(ny,1);
end

if(~exist('thresh', 'var') || isempty(thresh))
    thresh = 1e-3;
end

if(~exist('vectorized', 'var'))
    vectorized = 0;
end

if numel(thresh) == 1
   thresh = thresh(ones(size(y)));
end

 
yscale = max(abs(y),thresh);
del = (y + fac .* yscale) - y;
for j = find(del == 0)'
  while 1
    if fac(j) < facmax
      fac(j) = min(100*fac(j),facmax);
      del(j) = (y(j) + fac(j)*yscale(j)) - y(j);
      if del(j)
        break
      end
    else
      del(j) = thresh(j);
      break;
    end
  end
end
if nF == ny
  s = (sign(Fy) >= 0);
  del = (s - (~s)) .* abs(del);         % keep del pointing into region
end

% Form a difference approximation to all columns of dFdy.
ydel = y(:,ones(1,ny)) + diag(del);
if vectorized
    Fdel = feval(F,ydel);
else
    Fdel = zeros(nF,ny);
    for j = 1:ny
        Fdel(:,j) = feval(F,ydel(:,j));
    end
end
Fdiff = Fdel - Fy(:,ones(1,ny));
dFdy = Fdiff * diag(1 ./ del);
[Difmax,Rowmax] = max(abs(Fdiff),[],1);
% If Fdel is a column vector, then index is a scalar, so indexing is okay.
absFdelRm = abs(Fdel((0:ny-1)*nF + Rowmax)); 

% Adjust fac for next call to numjac.
absFy = abs(Fy);
absFyRm = absFy(Rowmax);              % not a col vec if absFty scalar
absFyRm = absFyRm(:)';                % ensure that absFtyRm is a row vector
absFdelRm = absFdelRm(:)';              % ensure that absFdelRm is a row vector
j = ((absFdelRm ~= 0) & (absFyRm ~= 0)) | (Difmax == 0);

if any(j)
  ydel = y;
  Fscale = max(absFdelRm,absFyRm);

  % If the difference in f values is so small that the column might be just
  % roundoff error, try a bigger increment. 
  k1 = (Difmax <= br*Fscale);           % Difmax and Fscale might be zero
  for k = find(j & k1)
    tmpfac = min(sqrt(fac(k)),facmax);
    del = (y(k) + tmpfac*yscale(k)) - y(k);
    if (tmpfac ~= fac(k)) && (del ~= 0)
      if nF == ny
        if Fy(k) >= 0                  % keep del pointing into region
          del = abs(del);
        else
          del = -abs(del);
        end
      end
        
      ydel(k) = y(k) + del;
      fdel = feval(F,ydel);
      ydel(k) = y(k);
      fdiff = fdel(:) - Fy;
      tmp = fdiff ./ del;
      
      [difmax,rowmax] = max(abs(fdiff));
      if tmpfac * norm(tmp,inf) >= norm(dFdy(:,k),inf)
        % The new difference is more significant, so
        % use the column computed with this increment.
        dFdy(:,k) = tmp;
  
        % Adjust fac for the next call to numjac.
        fscale = max(abs(fdel(rowmax)),absFy(rowmax));
          
        if difmax <= bl*fscale
          % The difference is small, so increase the increment.
          fac(k) = min(10*tmpfac, facmax);
          
        elseif difmax > bu*fscale
          % The difference is large, so reduce the increment.
          fac(k) = max(0.1*tmpfac, facmin);

        else
          fac(k) = tmpfac;
            
        end
      end
    end
  end
  
  % If the difference is small, increase the increment.
  k = find(j & ~k1 & (Difmax <= bl*Fscale));
  if ~isempty(k)
    fac(k) = min(10*fac(k), facmax);
  end

  % If the difference is large, reduce the increment.
  k = find(j & (Difmax > bu*Fscale));
  if ~isempty(k)
    fac(k) = max(0.1*fac(k), facmin);
  end
end
