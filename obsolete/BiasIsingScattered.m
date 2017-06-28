%run simulations for within-model class bias of Ising models:

if 0
    clear all
    addpath([pwd,'/../../FitMaxEnt']);
    addpath ~/CopyMPI/Code/minFunc/;
    basename='~/CopyMPI/MaxEnt/Data/BiasIsingScattered';
    thisdir=pwd;
    cd('~/CopyMPI/matlab/functions'), pathchoices=3; startup, cd(thisdir)

    mkdir(basename)
    callfile=mfilename;
    dims=[2,3,5,10,15]; %use Ising model with 2, 3 5, 10, 15 dimensions
    
    %means=[.05,.1,.5]; %try out means 0.01, .1, .5
    stdJs0=[.01,.1,1]; %standard devation of Js for population size 1, for population size n, std is divided by sqrt(n)
    numrethrowparams=10; %for each experiment, try 5 random throws of parameters J
    Ns=[10:10:200,round(10.^linspace(2,5,25))]; %number of 'datapoints' in simulation.
    Ns=unique(sort(Ns));
    
    numrepeats=10000;
    means=[.1,.5];
    
    
    
    %for each setting of parameters, create many data-sets to average over, so that we can calculate the expected bias, and the variance of the estimator
    numrepeats_save=5; %save detailed results for the first 10 data-sets:
    minmean=1e-20; %minimum allowed mean, any smaller mean is trunated at this value (otherwise maxent model is not defined, or lives in lower-dimensional space)
    
    %options for parameter learning:
    fitoptions.TolX=1e-50;
    fitoptions.TolFun=1e-50;
    fitoptions.display='off';
    fitoptions.MaxIter=5000;
    fitoptions.restarts=2;
    
    %first, set up all the parameters and save them:
    savenameparams=[basename,'Params'];
    
    for k=1:numel(means);
        for kk=1:numel(dims)
            for kkk=1:numel(stdJs0);
                for kkkk=1:numrethrowparams
            p.fmeanorig=means(k);
            p.dim=dims(kk);
            p.stdJs0=stdJs0(kkk);
            p.rethrow=kkkk;
            J=randn(p.dim)*p.stdJs0/sqrt(p.dim);
            h=log(p.fmeanorig/(1-p.fmeanorig));
            p.h=ones(p.dim,1)*h;
            p.J=sparse(J);
            p.lambda=hJ2lambda(p.h,p.J);
            [logP,p.logZ,P, p.means]=PMaxEnt(p.dim,p.lambda);
            p.EntropyLogE=ent(P);
            p.EntropyLog2=p.EntropyLogE*log2(exp(1));
            p.started=nan;
            p.finished=nan;
            p.filename=sprintf('%s/Results_%03d_%03d_%03d_%03d',basename,k,kk,kkk,kkkk);
            params(k,kk,kkk,kkkk)=p;
            %need checker that tells us when there are problems!!!!!
                end
            end
        end
    end
    clear log* P* h J
    save([basename,'/HQ']);
else
    %clear all
    thisdir=pwd;
    cd('~/CopyMPI/matlab/functions'), pathchoices=4; startup, cd(thisdir)
    addpath([pwd,'/../../FitMaxEnt']);
    addpath ~/CopyMPI/Code/minFunc/;
    basename='~/CopyMPI/MaxEnt/Data/BiasIsingScattered';
    load([basename,'/HQ']);
end
%%
%break



%do everything for first parameter setting, as a test case:
%starti=1;
%starti=10;
if ~exist('docleanup','var')
    docleanup=true;
end
if ~exist('starti','var')
    starti=1;
end
if ~exist('stoprep','var')
    stoprep=2000;
end


for i=starti:numel(params);
    
    tic
    p=params(i);
    p.started=now;
    results=struct;
    [features,description,x]=SetupFeaturesMaxEnt(p.dim,2);
    
    fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('\n Running parameter set %03d',i)
    fprintf('\nRunning simulations with parameters');
    fprintf('\ndimension: %02d, h: %02d stdJ0: %02d, rethrow: %02d\n',p.dim, round(100*mean(p.h))/100, p.stdJs0, p.rethrow)
    if exist([p.filename,'.mat'],'file')
        %try loading it
        A=load([p.filename]);
        results=A.results;
        
        if numel(results(end).ExitFlag)==numrepeats;
            fprintf('\n File already exists and finished, skipping...')
            continue
        elseif docleanup==false
            fprintf('\n File staretd, but too lazy to cleanup, moving on\n');
            continue
        else
            kkstart=numel(results(end).ExitFlag)+1;
            fprintf('\n File already exists, starting at repeat %d\n', kkstart)
        end
    else
        kkstart=1;
    end
    
    [logPtrue,logZtrue,Ptrue, meanstrue]=PMaxEnt(features,p.lambda,p.logZ);
    for kk=kkstart:stoprep
        if mod(kk,100)==1
            tic
            fprintf('Repeat %03d of %d  ', kk, numrepeats);
        end
        
        randn('seed',kk+100*k+1000*i); rand('seed',kk+100*k+1000*i);
        
        featuressampling_all=SampleDiscrete(features,Ptrue,max(Ns));
        
        for k=1:numel(Ns);
            N=Ns(k);
            %PrintStar(k+1)
            results(k).N=N;
            
            fitoptions.TolX=(1e-40);
            fitoptions.TolFun=(1e-40);
            fitoptions.display='off';
            fitoptions.MaxIter=5000;
            
            %tic
            %calculate true probability distribution:
            
            %sample from it:
            %calculate means of features from samples (waste of memory!!!)
            featuressampling=featuressampling_all(1:Ns(k),:);
            meanssampling=mean(featuressampling,1); clear featuressampled
            
            %truncate means to avoid zero-counts:
            truncator=zeros(1,numel(meanssampling));
            truncator(1:p.dim)=minmean;
            truncator(p.dim+1:end)=minmean^2;
            meanssamplingtruncated=max(meanssampling, truncator);
            meanssamplingtruncated=min(meanssamplingtruncated, 1-truncator);
            %fitoptions.lambda0=p.lambda;
            
            [lambdalearned,logZlearned, Plearned, meanslearned,output]=FitMaxEntLinear(features,meanssamplingtruncated, fitoptions);
            ExitFlag=output.exitflag;
            Accuracy=max(abs(meanslearned-meanssampling));
            ZerosInMean=sum(meanssampling==0);
            
            EntropyLogE=ent(exp(Plearned));
            EntropyLog2=EntropyLogE*log2(exp(1));
            
            if kk<=numrepeats_save
                details.EntropyLogE=EntropyLogE;
                details.EntropyLog2=EntropyLog2;
                details.lambda=lambdalearned;
                details.means=meanssampling;
                details.logZ=logZlearned;
                details.meanslearned=meanslearned;
                details.Accuracy=Accuracy;
                details.ExitFlag=ExitFlag;
                details.output=output;
                details.ZerosInMean=ZerosInMean;
                results(k).details(kk)=details;
                
            end
            %asdfsd
            
            results(k).EntropyLogE(kk)=EntropyLogE;
            results(k).EntropyLog2(kk)=EntropyLog2;
            results(k).ExitFlag(kk)=ExitFlag;
            results(k).Accuracy(kk)=Accuracy;
            results(k).ZerosInMean(kk)=ZerosInMean;
            results(k).fitoptions=fitoptions;
            
            %        keyboard
        end
        %fprintf('\n Done with this round, moving on to next N\n')
        if mod(kk,100)==0
            toc
            save(p.filename,'p','results','fitoptions')
            %fprintf('\n')
        end
        %toc
        
    end
    if kk==numrepeats
        p.finished=now;
        params(i)=p;
        % keyboard
        p.usetrueinitial=true;
        fprintf('\n Successfully saved results, next round!\n')
    else
        fprintf('\n Getting bored, save and move on to next round!\n')
    end
    save(p.filename,'p','results','fitoptions')
        
    %toc
    %catch
    %    warning('Something terrible happend, moving on...')
    %end
    
end

%% added some stuff


