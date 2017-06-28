%run simulations for within-model class bias of Ising models:


addpath([pwd,'/../../FitMaxEnt']);
addpath ~/CopyMPI/Code/minFunc/;
addpath ~/CopyMPI/Code/minFunc/;
addpath ~/CopyMPI/MaxEnt/matlab/MaxEntBias/functions
addpath ~/CopyMPI/BinCorr/matlab/FlatModels
addpath ~/CopyMPI/BinCorr/matlab/FlatDG
addpath ~/CopyMPI/BinCorr/matlab/functions/
addpath ~/CopyMPI/BinCorr/matlab/DGfitting/
thisdir=pwd;
cd('~/CopyMPI/matlab/functions'), pathchoices=4; startup, cd(thisdir)


if 0
    clear all
    basename='~/CopyMPI/MaxEnt/Data/ExploreMutualInformation_1';
    
    mkdir(basename)
    callfile=mfilename;
    dims=[5];
    stimdims=[1,5];
    meanstims=[-3,-2,-1];
    receptivefieldstds=[.1,1];
    noisevar=1;
    noisecorr=.1;
    T=1000;
    
    Ns=[10,20,30,40:20:200,300,500,1000,2500,5000]; %number of 'datapoints' in simulation.
    Ns=unique(sort(Ns));
    numrepeats=5000;
    
    numrepeats_save=3; %save detailed results for the first 3 data-sets:
    minmean=1e-20; %minimum allowed mean, any smaller mean is trunated at this value (otherwise maxent model is not defined, or lives in lower-dimensional space)
    
    %options for parameter learning:
    fitoptions.TolX=1e-50;
    fitoptions.TolFun=1e-50;
    fitoptions.display='off';
    fitoptions.MaxIter=5000;
    fitoptions.restarts=2;
    
    %first, set up all the parameters and save them:
    savenameparams=[basename,'Params'];
    
    for k=1:numel(dims);
        [features,description,x]=SetupFeaturesMaxEnt(dims(k),2);
        for kk=1:numel(stimdims);
            for kkk=1:numel(meanstims);
                for kkkk=1:numel(receptivefieldstds);
                    k
                    kk
                    kkk
                    kkkk
                    tic
                    p.dim=dims(k);
                    p.stimdim=stimdims(kk);
                    p.noisecorr=noisecorr;
                    p.meanstim=meanstims(kkk);
                    p.receptivefieldstd=receptivefieldstds(kkkk);
                    p.receptivefields=randn(p.dim,p.stimdim)*p.receptivefieldstd;
                    
                    p.stim=randn(p.stimdim,T);
                    p.input=p.receptivefields*p.stim+p.meanstim;
                    p.noisecovin=(eye(p.dim)*(1-p.noisecorr)+ones(p.dim)*p.noisecorr);
                    p.totalcovin=p.noisecovin+p.receptivefields*p.receptivefields';
                    p.Lambda=p.totalcovin;
                    p.gamma=ones(p.dim,1)*p.meanstim;
                    
                    [p.totalmeanbin, p.totalcorrbin]=Lambda2Rho(p.gamma, p.Lambda);
                    
                    for t=1:T
                        PrintStar(t)
                        [p.condmean(:,t),p.condcorr(t,:,:)]=Lambda2Rho(p.input(:,t),p.noisecovin);
                        p.P_DG(:,t)=ProbsDG(x,p.input(:,t),p.noisecovin);
                        p.Entropy_DG(t)=ent(p.P_DG(:,t));
                        
                        mean_features_DG=features'*p.P_DG(:,t);
                        [lambda_q2,P.logZ_q2(t), p.P_q2(:,t)]=FitMaxEntLinear(features,mean_features_DG', fitoptions);
                        p.P_q2(:,t)=exp(p.P_q2(:,t));
                        p.Entropy_q2(t)=ent(p.P_q2(:,t));
                        %fit maximum entropy model to each
                        %stimulus-condition response:
                    end
                    
                    p.P_DG_total=ProbsDG(x,p.gamma,p.Lambda);
                    p.Entropy_DG_total=ent(p.P_DG_total);
                    
                    p.P_DG_marginal=mean(p.P_DG,2);
                    p.Entropy_DG_marginal=ent(p.P_DG_marginal);
                    
                    
                    mean_features_DG=features'*p.P_DG_marginal;
                    [lambda_q2,P.logZ_q2_total, p.P_q2_total]=FitMaxEntLinear(features,mean_features_DG', fitoptions);
                    p.P_q2_total=exp(p.P_q2_total);
                    p.Entropy_q2_total=ent(p.P_q2_total);
                    
                    
                    p.P_q2_marginal=mean(p.P_q2,2);
                    p.Entropy_q2_marginal=ent(p.P_DG_marginal);
                    
                    p.MI_DG_total=p.Entropy_DG_total-mean(p.Entropy_DG);
                    p.MI_DG_marginal=p.Entropy_DG_marginal-mean(p.Entropy_DG);
                    p.MI_q2_total=p.Entropy_q2_total-mean(p.Entropy_q2);
                    p.MI_q2_marginal=p.Entropy_q2_marginal-mean(p.Entropy_q2);
                    
                    % keyboard
                    
                    
                    p.filename=sprintf('%s/Results_%03d_%03d_%03d_%03d',basename,k,kk,kkk,kkkk);
                    params(k,kk,kkk,kkkk)=p;
                    %need checker that tells us when there are problems!!!!!
                    toc
                end
                save([basename,'/HQ']);
            end
        end
    end
else
    %clear all
    basename='~/CopyMPI/MaxEnt/Data/ExploreMutualInformation_1';
    load([basename,'/HQ']);
    numrepeats_save=4;
    Ns=[10,20,30,40:20:200,300,500,1000,2500,5000]; %number of 'datapoints' in simulation.
    Ns=unique(sort(Ns));
    
end

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
    stoprep=1;
end
stoprep=5;

for i=starti:numel(params);
    
    tic
    p=params(i);
    p.started=now;
    results=struct;
    [features,description,x]=SetupFeaturesMaxEnt(p.dim,2);
    
    fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('\n Running parameter set %03d',i)
    fprintf('\nRunning simulations with parameters');
    fprintf('\ndim: %02d, stim dim: %02d mean stim: %02d, rec field std: %02d\n',p.dim, p.stimdim, p.meanstim, p.receptivefieldstd)
    if exist([p.filename,'.mat'],'file')
        %try loading it
        A=load([p.filename]);
        results=A.results;
        if numel(results(end).ExitFlag)==numrepeats;
            fprintf('\n File already exists and finished, skipping...')
            continue
        else
            kkstart=numel(results(end).ExitFlag)+1;
            fprintf('\n File already exists, starting at repeat %d\n', kkstart)
        end
    else
        kkstart=1;
    end
    
    for kk=kkstart:stoprep
        if mod(kk,10)==1
            tic
            fprintf('Repeat %03d of %d  ', kk, numrepeats);
        end
        
        randn('seed',kk+100*k+1000*i); rand('seed',kk+100*k+1000*i);
        
        for t=1:T
            featuressampling_all{t}=sparse(SampleDiscrete(features,p.P_DG(:,t),max(Ns)));
        end
        
        for k=1:numel(Ns);
            %fprintf('Running data-set size %g',N);
            N=Ns(k);
            fprintf('\nRunning data-set size %g',N);
            
            %PrintStar(k+1)
            results(k).N=N;
            
            fitoptions.TolX=(1e-40);
            fitoptions.TolFun=(1e-40);
            fitoptions.display='off';
            fitoptions.MaxIter=5000;
            
            
            %calculate entropy for each t separately
            for t=1:T
                %PrintStar(t)
                
                featuressampling=full(featuressampling_all{t}(1:Ns(k),:));
                meanssampling=mean(featuressampling,1); clear featuressampled
                truncator=zeros(1,numel(meanssampling));
                truncator(1:p.dim)=minmean;
                truncator(p.dim+1:end)=minmean^2;
                meanssamplingtruncated=max(meanssampling, truncator);
                meanssamplingtruncated=min(meanssamplingtruncated, 1-truncator);
                [lambdalearned,logZlearned, P_q2_learned(:,t), meanslearned,output]=FitMaxEntLinear(features,meanssamplingtruncated, fitoptions);
                P_q2_learned(:,t)=exp(P_q2_learned(:,t));
                
                Accuracy(t)=max(abs(meanslearned-meanssampling));
                ZerosInMean(t)=sum(meanssampling==0);
                EntropyLogE(t)=ent(P_q2_learned(:,t));
                EntropyLogECorrected(t)=EntropyLogE(t)+numel(lambdalearned)/2/N;
            end
            
            P_q2_marginal=mean(P_q2_learned,2);
            Entropy_q2_marginal=ent(P_q2_marginal);
            
            
            meanssampling=(features'*P_q2_marginal)';
            truncator=zeros(1,numel(meanssampling));
            truncator(1:p.dim)=minmean;
            truncator(p.dim+1:end)=minmean^2;
            meanssamplingtruncated=max(meanssampling, truncator);
            meanssamplingtruncated=min(meanssamplingtruncated, 1-truncator);
            
            fitoptions.lambda0=lambdalearned;
            [lambda_q2,logZ_q2_total, P_q2_total,fitmeans,output]=FitMaxEntLinear(features,meanssamplingtruncated, fitoptions);
            P_q2_total=exp(P_q2_total);
            Entropy_q2_total=ent(P_q2_total);
            %keyboard
            
            results(k).EntropyLogE(kk,:)=EntropyLogE;
            results(k).CondEntropyLogE(kk,1)=mean(EntropyLogE(:));
            results(k).EntropyLogECorrected(kk,:)=EntropyLogECorrected;
            results(k).CondEntropyLogECorrected(kk,1)=mean(EntropyLogECorrected(:));
            
            
            
            results(k).Entropy_q2_total(kk,1)=Entropy_q2_total;
            results(k).Entropy_q2_total_Corrected(kk,1)=Entropy_q2_total+numel(lambdalearned)/2/N/T;
            
            results(k).Entropy_q2_marginal(kk,1)=Entropy_q2_marginal;
            results(k).Entropy_q2_marginal_Corrected(kk,1)=Entropy_q2_marginal+numel(lambdalearned)/2/N/T;
            
            results(k).MI_q2_total(kk)=results(k).Entropy_q2_total(kk)-results(k).CondEntropyLogE(kk);
            results(k).MI_q2_marginal(kk)=results(k).Entropy_q2_marginal(kk)-results(k).CondEntropyLogE(kk);
            
            results(k).MI_q2_total_Corrected(kk)=results(k).Entropy_q2_total_Corrected(kk)-results(k).CondEntropyLogECorrected(kk);
            results(k).MI_q2_marginal_Corrected(kk)=results(k).Entropy_q2_marginal_Corrected(kk)-results(k).CondEntropyLogECorrected(kk);
            
            results(k).ExitFlag(kk)=1;
            % results(k)
            %keyboard
            
        end
        if mod(kk,1)==0
            toc
            save(p.filename,'p','results','fitoptions')
        end
        
    end
    results
    if kk==numrepeats
        p.finished=now;
        params(i)=p;
        p.usetrueinitial=true;
        fprintf('\n Successfully saved results, next round!\n')
    else
        fprintf('\n Getting bored, save and move on to next round!\n')
    end
    save(p.filename,'p','results','fitoptions')
    
end




