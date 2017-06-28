
function [xSampled] = maxEnt_gibbs_fat(nSamples, burnIn, thinning, lambda, x0, model)
% This is the old version, called 'fat' to distinguish it from the 'slim'
% version that's trimmed to be faster, but only on 'ising_count_l_0' models
% input:
% -nSamples: number of Gibbs samples to be generated
% -  burnIn: number of to-be-discarded samples from beginning of MCMC chain
% -thinning: number of samples to be discarded between each stored pair
% -  lambda: parameter vector for maxEnt model
% -      x0: EITHER a d-by-1 vector for initial member of the MCMC,
%            OR a single number specifing the dimensionality of data d
%            In the latter case, initial chain member will be drawn.
% -   model: string specifying the layout of the feature function for the
%            maxEnt model. Supported values are "ising" (pure Ising model),
%            "ising_count" with additional terms for activity counts V(K),
%            "ising_count_l_0" with yet an adiditonal term for V(K) = 0.

if numel(x0) == 1
 d = x0;                                      % Generate x0 using E[X] from
 EX = exp(lambda(1:d))./(1+exp(lambda(1:d))); % a maxEnt model with only h,
 x0 = double(rand(d,1)<EX);                   % i.e. no parameters J, L
end
d = length(x0);

xSampled = zeros(d, nSamples); % +1 for x0

fxc = maxEnt_features_slim(x0, model);  % mostly needed to get dimensionality
m   = zeros(length(fxc),d); 
for k = 1:d % compute masks
   m(k,k)= 1;           % add terms corresponding to features for h   
   tmp = zeros(d,d);    % add terms corresponding to features for J
   tmp(k,:) = 1; tmp(:,k) = 1;
   tmp = tmp(logical(tril(ones(d,d),-1)));
   m(d+(1:d*(d-1)/2),k)   = tmp(:);        
end

xc = x0; % current sample
for i = 1:thinning*nSamples+burnIn 
    
 k = randi(d); % might choose another template for going through dimensions
                 
  % Compute probability
   v = m(:,k); % current mask for choosing correct feature components. 
               % Last entries depend on current sample xc for 
               % model = 'ising_count', 'ising_count_l_0'
   % compute p(x_k = 0), finisch mask in case of terms for V(K) are needed
   switch model
       case 'ising'
        p0 = 1;              
       case 'ising_count'
        idl = sum(xc) - xc(k);           
        v(d*(d+1)/2+(idl+1)) = 1;
        if idl>0, p0 = exp(-lambda(end-d+1+idl)); else p0 = 1; end
       case 'ising_count_l_0'
        idl = sum(xc) - xc(k);           
        v(d*(d+1)/2+(idl+2)) = 1;
        p0 = exp(-lambda(end-d+idl));
   end
   % compute p(x_k = 1)
   xc1 = xc; xc1(k) = 1;
   fxc1 = maxEnt_features_slim(xc1, model); % features, needed for v for V(K)
   p1 = exp(-(lambda.*v)' * fxc1);
   p1  = p1 / (p0 + p1); % normalization step
 
  
  % Update chain  
    xc(k) = double(rand(1) < p1);  % Draw new k-th entry
    if i>burnIn && mod(i-burnIn, thinning)==0
     xSampled(:,(i-burnIn)/thinning) = xc;
    end
end
 
end

function fx = maxEnt_features_slim(x, model)
 d = length(x);
 switch model
    case {2,'ising'}
        pairs=nchoosek(1:d,2);
        %description=[[1:d; (1:d)*nan],pairs];  
        fx=[x;x(pairs(:,1)).*x(pairs(:,2))];
    case 3
        pairs=nchoosek(1:d,2);
        triplets=nchoosek(1:d,3);
        %description=[[1:d; (1:d)*nan;(1:d)*nan],[pairs; pairs(1,:)*nan],triplets];  
        fx=[x;x(pairs(:,1)).*x(pairs(:,2));x(pairs(:,1)).*x(pairs(:,2)).*x(triplets(:,3))];
    case 'ising_count'
        pairs=nchoosek(1:d,2);
        count_indicators=zeros(d,1);
        sum_x = sum(x);
        if sum_x>1
         count_indicators(sum_x)=1;
        end
        %description=[[1:d; (1:d)*nan],pairs,[1:d; (1:d)*nan]];  
        fx=[x;x(pairs(:,1)).*x(pairs(:,2));count_indicators];
        %keyboard
    case 'ising_count_l_0' % same as 'ising_count', but has feature for K=0
        pairs=nchoosek(1:d,2);
        count_indicators=zeros(d+1,1);
        count_indicators(sum(x)+1) = 1; 
        %description=[[1:d; (1:d)*nan],pairs,[0:d; (0:d)*nan]];  
        fx=[x;x(pairs(:,1)).*x(pairs(:,2));count_indicators];
end

end
