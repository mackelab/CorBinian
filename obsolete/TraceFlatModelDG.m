
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
basename='~/CopyMPI/MaxEnt/Data/TraceFlatModelDG';
    

        mkdir(basename)
    %end
    callfile=mfilename;
    
    load([basename,'/allstuff']);
    
end


dims=[4,5,6,8,10:2:50,55:5:200]; %use Ising model with 2, 5, 10, 15 dimensions
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
            
            everythingok=  exist('params','var') && k<=size(params,1) && kk<=size(params,2) && kkk<=size(params,3) && (~isempty(params(k,kk,kkk)) && numel(params(k,kk,kkk).model)==5);
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
            p.model{3}=SampleFromFlatModel('maxbias',p.meanorig,p.corrorig,p.dim,1,1,ops);
            [p.model{3}]=CalcBxEtcFlatModel(p.model{3}, p.model{1});
            p.model{3}.S=[];
            
            ops.Bias=p.model{1}.Bias*1.01;
            p.model{4}=SampleFromFlatModel('maxbias',p.meanorig,p.corrorig,p.dim,1,1,ops);
            [p.model{4}]=CalcBxEtcFlatModel(p.model{4}, p.model{1});
            p.model{4}.S=[];
            
            ops.Bias=p.model{1}.Bias*1.1;
            p.model{5}=SampleFromFlatModel('maxbias',p.meanorig,p.corrorig,p.dim,1,1,ops);
            [p.model{5}]=CalcBxEtcFlatModel(p.model{5}, p.model{1});
            p.model{5}.S=[];
            
            
            % keyboard
            for kkkk=1:5
                p.Entropy(kkkk)=p.model{kkkk}.H_ana;
                p.Bias(kkkk)=p.model{kkkk}.Bias;
                % p.model{kkkk}.Sigma=p.model{kkkk}.S.Sigma;
                p.model{kkkk}.S=[];
                p.model{kkkk}.model_q2=[];
                p.model{kkkk}.p=[];
                p.model{kkkk}.Bx_canonical=single(p.model{kkkk}.Bx_canonical);
                p.FitError(kkkk)=p.model{kkkk}.error;
            end
            p.DeltaS=[0,p.Entropy(1)-p.Entropy(2),p.Entropy(1)-p.Entropy(3),p.Entropy(1)-p.Entropy(4),p.Entropy(1)-p.Entropy(5)];
            p.bDash=p.model{1}.bDash;
            
            p.DeltaS(6)=(p.Bias(3)-p.Bias(1)).^2/p.bDash/2;
            p.Entropy(6)=p.Entropy(1)-p.DeltaS(6);
            
            p.DeltaS(7)=(p.Bias(4)-p.Bias(1)).^2/p.bDash/2;
            p.Entropy(7)=p.Entropy(1)-p.DeltaS(7);
            
            p.DeltaS(8)=(p.Bias(5)-p.Bias(1)).^2/p.bDash/2;
            p.Entropy(8)=p.Entropy(1)-p.DeltaS(8);
            
            %if min(p.DeltaS(4:end))<-1e-5;
            %    keyboard
            %end
            
            params(k,kk,kkk)=p;
            toc
            %clear junk
        end
        if k>30
                   save([basename,'/allstuff'])
        end
    end
    if k>27
    save([basename,'/allstuff'])
    end
end
%%
%at the moment:
%find maximal

%plots that I want, as function of parameters mu and rho:
%
%additional bias of DG
%additional bias of maximally biased model
%perturbative prediction of maximally biased model
%bdash as function of parameters







