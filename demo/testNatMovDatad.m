clear all;
close all
load(['C:\Users\Loki\Desktop\Criticality\Results\run_40_735904.mat']);
%load('/home/marcel/criticality/results/NatMovData.mat')
nComp=40;
pairs=nchoosek(1:nComp,2);
nSamplesEval = 25000;
burnIn = 5000;
fracs = [0.5,2,Inf];
clrs = hsv(length(fracs));
j = 1;
lgnd = cell(length(fracs),1);
x0 = randi(2,[nComp,1])-1;
for i = 1:length(fracs)
      disp([num2str(i), ' of ' , num2str(length(fracs))])
        lambdaMPF = results.lambda;
%        lambdaMPF(nComp+1:nComp*(nComp+1)/2)=lambdaMPF(nComp+1:nComp*(nComp+1)/2)/fracs(i);
        lambdaMPF(end-nComp:end)=lambdaMPF(end-nComp:end)/fracs(i);
        [mfxEval,~,xc] = maxEnt_gibbs_pair_C(nSamplesEval, burnIn, lambdaMPF, x0, 'laptop');
        xc'
        figure(1)
        subplot(1,2,1)
         plot(lambdaMPF, 'color', clrs(i,:)); hold on;
        subplot(1,2,2)
         plot(mfxEval, 'color', clrs(i,:)); hold on;
        lgnd{i} = ['fraction ', num2str(fracs(i))];
        
        Ceval = zeros(nComp,nComp);
        for j= 1:size(pairs,1)
            k = pairs(j,1); l = pairs(j,2);
           Ceval(k,l) = mfxEval(j+nComp)-mfxEval(k)*mfxEval(l);
        end
        for j = 1:nComp
            Ceval(j,j) = mfxEval(j)*(1-mfxEval(j));
        end
        Ceval = Ceval + Ceval' - diag(diag(Ceval));
        CorEval = zeros(nComp,nComp);
        for k = 1:nComp
            for l=(k+1):nComp
                CorEval(k,l) = Ceval(k,l)/(sqrt(Ceval(k,k))*sqrt(Ceval(l,l)));
            end
        end
        CorEval = CorEval + CorEval' - diag(diag(CorEval));
        figure(2)
        subplot(2,3,i),
        title(['fraction ', num2str(fracs(i))])
        imagesc(Ceval); colorbar
clear k l i

end
figure(1)
subplot(1,2,1)
plot(results.lambda,'k.-')
subplot(1,2,2)
plot(results.mfxTrain, 'k.-')
disp('done')
legend(lgnd);
figure(2)
subplot(2,3,6)
imagesc(data.spkCov); colorbar
title('true')
%imagesc(data.spkCorrs-diag(diag(data.spkCorrs))); colorbar