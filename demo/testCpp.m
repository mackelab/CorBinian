nSamples = 5000; 
burnIn=100; 
d = 25; 
xc = randi(2,[d,1])-1; 
pairs = nchoosek(1:d,2); 
h = randn(d,1);  
J = randn(d*(d-1)/2,1)/sqrt(d); 
L = randn(d+1,1); 

m   = false(d*(d-1)/2,d); % m indexes in d-by-d matrices which d-1 entries
fm =  zeros(d-1, d);      % on the upper diagonal half share a particular
% index k, in booleans. fm = find(m), per column.
for k = 1:d
tmp = false(d,d);
tmp(k,:) = true; tmp(:,k) = true;
tmp = tmp(logical(tril(ones(d,d),-1)));
m(:,k)   = tmp(:);
fm(:,k) = find(m(:,k));
end
clear tmp
%%
tic 
[xSampledC] = pwGibbsMaxEnt_malloc(int32(nSamples), int32(burnIn), int32(d), ...
                    xc, pairs-1, m, fm-1, h, J, L);
toc
%%
tic
[xSampled, E] = GibbsLoop(nSamples,burnIn,d,xc,pairs,m,fm,h,J,L);
toc
%%
%[xSampledC,xSampled];
%%
figure(3); 
plot(xSampled, 'go-'), hold on; plot(xSampledC, 'ro-'), hold off;