% this script serves to test the equations derived for blocked iterative
% scaling on the full maximum entropy model for neural spike patterns
% devised by W. Bialek & co.
newLambda = true;
newData   = true; nSamplesTrain = 100000; burnIn = 10000; 
n = 10;

if newLambda
 h=randn(n,1)-1;  J=randn(n); J=triu(J,1)/sqrt(n); lambda=hJ2lambda(h,J);
 L=randn(n+1,1)/sqrt(n); 
 L = L - L(1); % a rather recently discovered nasty extra degree of freedom
 lambdaTrue = [lambda;L];
 lambdaHat  = lambdaTrue + randn(size(lambdaTrue)); 
end

if newData
 disp('Generating training data')
 [Efx,~,~] = maxEnt_gibbs_pair_C(nSamplesTrain, burnIn, lambdaTrue, n);
 disp('Generating MCMC sample')
 [Efy,~,~] = maxEnt_gibbs_pair_C(nSamplesTrain, burnIn, lambdaHat,  n);
end
%%

for r = 1:100
idxij = nchoosek(1:n,2);
i = randsample(n,2);
j = max(i); i = min(i); ij = find((idxij(:,1)==i)&(idxij(:,2)==j)) + n; 

di = zeros(11,1); dj = zeros(size(di));
di(1) = 1; dj(1) = 1;
for t = 1:10
    di(t+1) = log(Efx(i)/(1-Efx(i))) + ...
         log( (1+(exp(dj(t))-1)*Efy(j))/(Efy(i)+(exp(dj(t))-1)*Efy(ij)) - 1);
    dj(t+1) = log(Efx(j)/(1-Efx(j))) + ...
         log( (1+(exp(di(t))-1)*Efy(i))/(Efy(j)+(exp(di(t))-1)*Efy(ij)) - 1);
end

%     di2 = log(Efx(i)/(1-Efx(i)))  ...
%         + log( (1-Efy(i)-Efy(j)+Efy(ij)) * ( (exp(di(end))-1)*Efy(ij)+Efy(j) ) + exp(di(end)) * Efx(j) *(Efy(i)*Efy(j)-Efy(ij))) ...
%         - log( (   Efy(i) - Efy(ij)    ) * ( (exp(di(end))-1)*Efy(ij)+Efy(j) ) -                Efx(j) *(Efy(i)*Efy(j)-Efy(ij)))
%     dj2 = log(Efx(j)/(1-Efx(j)))  ...
%         + log( (1-Efy(j)-Efy(i)+Efy(ij)) * ( (exp(dj(end))-1)*Efy(ij)+Efy(i) ) + exp(dj(end)) * Efx(i) *(Efy(i)*Efy(j)-Efy(ij))) ...
%         - log( (   Efy(j) - Efy(ij)    ) * ( (exp(dj(end))-1)*Efy(ij)+Efy(i) ) -                Efx(i) *(Efy(i)*Efy(j)-Efy(ij)))


 a = (1-Efx(i)) * Efy(ij) * (Efy(i)-Efy(ij));
 b = (Efx(i)+Efx(j))*(Efy(ij)-Efy(i)*Efy(j))+(Efy(ij)-Efy(i))*(Efy(ij)-Efy(j));
 c = Efx(i)*(Efy(j)-Efy(ij))*(1+Efy(ij)-Efy(i)-Efy(j));
 x = exp(di(end));
di2 = log( (a*x^2 + ( b + (Efx(j)*(Efy(i)*Efy(j)-Efy(ij))) ) * x + c) / (Efx(j)*(Efy(i)*Efy(j)-Efy(ij))) );

    plot(di); hold on; plot(dj, 'r')
    plot(1.2*t, di2, 'b*');
%    plot(1.2*t, dj2, 'r*');
    pause;
    hold off
end


% a(r) = (1-Efx(j)) * Efy(ij) * (Efy(j)-Efy(ij));
% b(r) = (Efx(i)+Efx(j))*(Efy(ij)-Efy(i)*Efy(j))+(Efy(ij)-Efy(i))*(Efy(ij)-Efy(j));
% c(r) = Efx(j)*(Efy(i)-Efy(ij))*(1+Efy(ij)-Efy(i)-Efy(j));
% 
% [a,b,c]
% 
% dia1(r) = (-b(r) - sqrt(b(r)^2 - 4*a(r)*c(r)))/(2*a(r));
% dia2(r) = (-b(r) + sqrt(b(r)^2 - 4*a(r)*c(r)))/(2*a(r));
% 
% %figure; 
%         plot(di); hold on; plot(dj, 'r')
%         plot(1.3*t, dia1(r), 'b*');
%         plot(1.3*t, dia2(r), 'bo');        
% 
% pause;
% end