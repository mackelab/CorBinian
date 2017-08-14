function [xSampled, E] = GibbsLoop(nSamples,burnIn,d,xc,pairs,m,fm,h,J,L)

xSampled = zeros(d*(d+3)/2+1, 1);      % return Bernoulli probabilities
xTmp = xSampled; % intermediate storage for normalization (numerical issue)
E = zeros(nSamples,1); % return sequence of energies
idl = sum(xc)+1;

for i = -burnIn:nSamples
 ks = randperm(d*(d-1)/2); % one MCMC update equals one sweep through all d data- 
                           % dimensions in random order
 pc = zeros(d*(d+3)/2+1,1);  % current Bernoulli probabilities
 for j = 1:length(ks)          
   k = pairs(ks(j),1); 
   l = pairs(ks(j),2);
   
   idl = idl - xc(k) - xc(l); % current activity count, x(k), xc(l) IGNORED                
  
  % compute some key values
   xc1 = xc; 
   xc1(k) = true; xc1(l) = false; 
   sk = sum( J( fm(xc1(pairs(m(:,k),1))&xc1(pairs(m(:,k),2)),k) ) );
   % m(:,k) holds which of the d-1 indices of J are involved with x(k)
   % pairs(m(:,k),1:2) gives the d-1 index pairs (i,j) of x that correspond
   %                   to the indices of J that are involved with x(k)
   % xc1(pairs(...,1))&x1(pairs(...,2)) filters out those where x(i)=x(j)=1
   % fm(xc1(...)&xc1(...)) translates the pair indices back into J indices
   
   xc1(k) = false; xc1(l) = true; 
   sl = sum( J( fm(xc1(pairs(m(:,l),1))&xc1(pairs(m(:,l),2)),l) ) );
   
   xc1(k) = true;  
   slk = sum( J( fm(xc1(pairs(m(:,k),1))&xc1(pairs(m(:,k),2)),k) ) );
   
   L1 = L(idl+1); % one of the two variables xc(k), xc(l) is = 1
   L2 = L(idl+2); % both variables are = 1
      
  % compute joint probabilities for bivariate Bernoulli   
   p = exp([h(k) +        sk  +       L1           % xc(k) = 1, xc(l) = 0
            h(k) + h(l) + sl  + slk + L2;          % xc(k) = 1, xc(l) = 1
                   h(l) +        sl + L1;          % xc(k) = 0, xc(l) = 1
                                      L(idl)]);    % xc(k) = 0, xc(l) = 0
   p = p/sum(p); % normalize
  % update chain   
   rnd = rand(1);               % essentially, use the CDF for each of the
   xc(k) = (rnd < p(1)+p(2));   % four possible outcomes for (xc(k),xc(l))
   xc(l) = (p(1) < rnd && rnd < 1-p(4));       
   idl = idl + xc(k) + xc(l);
   
   pc(k) = pc(k) + p(1) + p(2);  % for E[x(k)]
   pc(l) = pc(l) + p(2) + p(3);  % for E[x(l)]
   pc(d+ks(j)) = p(2); % for E[x(k) * x(l)] and thus for Cov(x(k),x(l))
   pc(d*(d+1)/2+idl) = pc(d*(d+1)/2+idl)+1; % for E[sum(x)==k], i.e. V(K)
 end % end for j=1:length(ks), i.e. loop over all updates within one sweep
 
 pc(1:d) = pc(1:d) / (d-1);  % we just visited each x(k) d-1 times
 pc(d*(d+1)/2+(1:d+1)) = pc(d*(d+1)/2+(1:d+1))/ length(ks);
 
 %E(i) = h'*xc + xc'*J*xc + L(idl);
 
 % Update online estimate of expected values
  if i>0
   xTmp = xTmp + pc; 
   if ~mod(i,1000)% normalize output every 1000 samples to avoid numerical 
    xSampled = ((i-1000)*xSampled + xTmp)/i; % issues with large numbers.
    xTmp = 0;
    disp(['  - ', num2str(i),'/',num2str(nSamples),'samples'])
   end
 end
      
end % end for i=1:nSamples

end