
%% %% %%

%simulate from Ising model and DG for increasing population size, and fit
%bias analytically and numerically in both cases, to see how the bias
%scales with population size:
%changed to 200

if 1
    clear all
    thisdir=pwd;
    cd('~/CopyMPI/matlab/functions'), pathchoices=4; startup, cd(thisdir);
    
    addpath([pwd,'/../../FitMaxEnt']);
    addpath ~/CopyMPI/Code/minFunc/;
    addpath ~/CopyMPI/MaxEnt/matlab/MaxEntBias/functions
    addpath ~/CopyMPI/BinCorr/matlab/FlatModels
    addpath ~/CopyMPI/BinCorr/matlab/FlatDG
    addpath ~/CopyMPI/BinCorr/matlab/functions/
    addpath ~/CopyMPI/BinCorr/matlab/DGfitting/
    % try
    basename='~/CopyMPI/MaxEnt/Data/MaxMinBiasModel';
    
    
    mkdir(basename)
    %end
    callfile=mfilename;
    
    load([basename,'/allstuff']);
    
end

scalers=linspace(0,1,11);
scalers=[scalers(2:end-1),exp(linspace(0,log(10),41))];
scalers=[- fliplr(scalers), 0, scalers];

%calculate minDeltaS model using Bias deviations which are
%+- 10 times the bias of the DG, in bins which increase in sqare steps

dims=[5,10,50,100]; %use Ising model with 2, 5, 10, 15 dimensions
means=[.02,.05,.1,.2,.5]; %try out means 0.01, .1, .5
rhos=[.01,.02,.05,.1,.2,.5]; %take pairwise correlation 0, .1, .5

%dims=[10,20]; %use Ising model with 2, 5, 10, 15 dimensions
%means=[.1:.1:.5]; %try out means 0.01, .1, .5
%rhos=[.1:.1:.5]; %take pairwise correlation 0, .1, .5





%first, set up all the parameters and save them:
%    savenameparams=[basename,'Params'];
randn('seed',0); rand('seed',0);
global mutarget vartarget;

for k=1:numel(dims);
    k
    N=dims(k)
    pause(1)
    for kk=1:numel(means);
        kk
        for kkk=1:numel(rhos);
            % try
            %for kkkk=1:5
            %    pspike=means(kk);
            %    rho=rhos(kkk);
            %    params(k,kk,kkk).model{kkkk}.mucount_target=N*pspike;
            %    params(k,kk,kkk).model{kkkk}.varcount_target=pspike*(1-pspike)*N*(1+(N-1)*rho);
            %
            %end
            %mutarget=N*means(kk);
            %vartarget=means(kk)*(1-means(kk))*N*(1+(N-1)*rhos(kkk));
            %keyboard
            everythingok=  exist('params','var') && k<=size(params,1) && kk<=size(params,2) && kkk<=size(params,3) && (~isempty(params(k,kk,kkk)) && numel(params(k,kk,kkk).model)==4);
            %if everythingok
            %    accurate=true;
            %    for kkkk=[1,3,4,5];
            
            %        acc=max(abs([mutarget-params(k,kk,kkk).model{kkkk}.mucount_ana,vartarget-params(k,kk,kkk).model{kkkk}.varcount_ana]));
            %        accurate=accurate && acc<1e-5;
            %  if ~accurate
            %      keyboard
            %  end
            %    end
            %end
            if everythingok %%&& accurate
                continue
                
            end
            % end
            tic
            %kkk
            p=struct;
            p.dim=dims(k);
            p.meanorig=means(kk);
            p.corrorig=rhos(kkk);
            
            p.model{1}=SampleFromFlatModel('ising',p.meanorig,p.corrorig,p.dim,1,1);
            [p.model{1}]=CalcBxEtcFlatModel(p.model{1}, p.model{1});
            p.model{1}.S=[];
            
            p.model{2}=SampleFromFlatModel('DG',p.meanorig,p.corrorig,p.dim,1,1);
            [p.model{2}]=CalcBxEtcFlatModel(p.model{2}, p.model{1});
            p.model{2}.S=[];
            
            ops.Bias=p.model{2}.Bias;
            p.model{3}=SampleFromFlatModel('minDeltaS',p.meanorig,p.corrorig,p.dim,1,1,ops);
            [p.model{3}]=CalcBxEtcFlatModel(p.model{3}, p.model{1});
            p.model{3}.S=[];
            p.model{3}.model_q2=[];
            
            %ops.DeltaS=(p.model{1}.H_ana-p.model{2}.H_ana)/log2(exp(1));
            %warning off
            %p.model{4}=SampleFromFlatModel('maxBias',p.meanorig,p.corrorig,p.dim,1,1,ops);
            %warning on
            %[p.model{4}]=CalcBxEtcFlatModel(p.model{4}, p.model{1});
            
            p.extraBiasDG=p.model{2}.Bias-p.model{1}.Bias;
            p.allBiases.Biases=p.model{1}.Bias+scalers*max(1e-6,p.extraBiasDG);
            p.allBiases.Biases=p.allBiases.Biases(p.allBiases.Biases>0);
            %keyboard
            
            
            % keyboard
            for kkkk=1:3
                p.Entropy(kkkk)=p.model{kkkk}.H_ana/log2(exp(1));
                p.Bias(kkkk)=p.model{kkkk}.Bias;
                % p.model{kkkk}.Sigma=p.model{kkkk}.S.Sigma;
                p.model{kkkk}.S=[];
                p.model{kkkk}.model_q2=[];
                p.model{kkkk}.p=[];
                p.model{kkkk}.Bx_canonical=single(p.model{kkkk}.Bx_canonical);
                p.model{kkkk}=RmEmptyFields(p.model{kkkk});
                p.FitError(kkkk)=p.model{kkkk}.error;
            end
            p.DeltaS=[0,p.Entropy(1)-p.Entropy(2),p.Entropy(1)-p.Entropy(3)];
            p.bDash=p.model{1}.bDash;
            
            p.DeltaS(4)=(p.Bias(3)-p.Bias(1)).^2/p.bDash/2;
            p.Entropy(4)=p.Entropy(1)-p.DeltaS(4);
            
            toc
            tic
            % keyboard
            warning off
            for ki=1:numel(p.allBiases.Biases)
                
                PrintStar(ki)
                try
                    ops.Bias=p.allBiases.Biases(ki);
                    locmodel=SampleFromFlatModel('minDeltaS',p.meanorig,p.corrorig,p.dim,1,1,ops);
                    locmodel.model_q2=[];
                    locmodel=RmEmptyFields(locmodel);
                    p.allBiases.FitError(ki)=locmodel.error;
                    p.allBiases.FitErrorBias(ki)=locmodel.Bias_error;
                    
                    p.allBiases.Bias(ki)=locmodel.Bias;
                    p.allBiases.DeltaS(ki)=p.model{1}.H_ana-locmodel.H_ana;
                catch
                    p.allBiases.FitError(ki)=nan;
                    p.allBiases.FitErrorBias(ki)=nan;
                    
                    p.allBiases.Bias(ki)=nan;
                    p.allBiases.DeltaS(ki)=nan;
                end
            end
            warning on
            %keyboard
            params(k,kk,kkk)=p;
            toc
            %clear junk
            clear p
        end
        save([basename,'/allstuff'])
    end
end
