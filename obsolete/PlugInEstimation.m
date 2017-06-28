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
    basename='~/CopyMPI/MaxEnt/Data/PlugInEstimation_1';
    
    mkdir(basename)
    %load([basename,'/HQ']);
    callfile=mfilename;
    dims=[3,5,8,10]; %use Ising model with 2, 3 5, 10, 15 dimensions
    means=[.05,.1,.5];
    rhos=[0,.05,.1,.2,.5];
    scatterstds=[0,.1,.5,1];
    
    Ns=[10:10:200,round(10.^linspace(2,5,25))]; %number of 'datapoints' in simulation.
    Ns=unique(sort(Ns));
    
    numrepeats=5000;
    
    
    
    
    %for each setting of parameters, create many data-sets to average over, so that we can calculate the expected bias, and the variance of the estimator
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
        for kk=1:numel(means);
            for kkk=1:numel(rhos);
                for kkkk=1:numel(scatterstds);
                    k
                    kk
                    kkk
                    kkkk
                    tic
                    p.dim=dims(k);
                    p.fmeanorig=means(kk);
                    p.startcorr=rhos(kkk);
                    p.scatterstd=scatterstds(kkkk);
                    %first, make random gamma vector and correlation matrix:
                    
                    Lambda0=(eye(p.dim)*(1-p.startcorr)+ones(p.dim)*p.startcorr);
                    G0=chol(Lambda0);
                    G=G0+triu(randn(size(G0)),1)*p.scatterstd/sqrt(p.dim);
                    p.Lambda=Cov2Corr(G'*G);
                    gamma0=norminv(ones(p.dim,1)*p.fmeanorig);
                    
                    p.gamma=gamma0+randn(p.dim,1)*p.scatterstd;
                    [p.mu, p.Rho]=Lambda2Rho(p.gamma, p.Lambda);
                    p.meanmean=mean(p.mu);
                    p.medmean=median(p.mu);
                    p.harmman=exp(mean(log(p.mu)));
                    p.stdmean=std(p.mu);
                    p.rangemean=max(p.mu)-min(p.mu);
                    
                    corrs=p.Rho(triu(ones(size(p.Rho)),1)==1);
                    abscorrs=abs(corrs);
                    
                    p.meancorr=mean(corrs);
                    p.medcorr=median(corrs);
                    p.harmcorr=exp(mean(log(corrs)));
                    p.stdcorr=std(corrs);
                    p.rangecorr=max(corrs)-min(corrs);
                    
                    p.meanabscorr=mean(abscorrs);
                    p.medabscorr=median(abscorrs);
                    p.harmabscorr=exp(mean(log(abscorrs)));
                    p.stdabscorr=std(abscorrs);
                    p.rangeabscorr=max(abscorrs)-min(abscorrs);
                    
                    
                    %DG first:
                    p.P_DG=ProbsDG(x,p.gamma,p.Lambda)';
                    mean_features_DG=features'*p.P_DG;
                    p.Entropy_DG=ent(p.P_DG);
                    
                    %then Ising model
                    [p.lambda_q2,P.logZ_q2, p.P_q2, p.mean_features_q2,output]=FitMaxEntLinear(features,mean_features_DG', fitoptions);
                    p.P_q2=exp(p.P_q2);
                    p.Entropy_q2=ent(p.P_q2);
                    
                    
                    %then Bias calculations
                    p.Sigma_q2=wcov(features,p.P_q2);
                    p.Bx=CalcBx(features, p.mean_features_q2,p.Sigma_q2);
                    p.Bias_q2=p.Bx'*p.P_q2;
                    p.Bias_DG=p.Bx'*p.P_DG;
                    featuresBiased=[features,p.Bx];
                    [p.bDash,Hx,Bx,VarB,EBdeltag,SigmaBdeltag]=CalcbDashAtZero(features,p.mean_features_q2, p.Sigma_q2, p.P_q2);
                    
                    %then MaxBias model
                    meansBiased=[p.mean_features_q2,p.Bias_DG];
                    try
                        [p.lambda_biased,p.logZ_biased, p.P_biased, p.means_biased,output]=FitMaxEntLinear(featuresBiased,meansBiased, fitoptions);
                        p.P_biased=exp(p.P_biased);
                        p.Accuracy=max(abs(p.means_biased-meansBiased));
                        p.Entropy_biased=ent(P_biased);
                    catch
                        % keyboard
                        p.lambda_biased=nan;
                        p.logZ_biased=nan;
                        p.P_biased=nan;
                        p.means_biased=nan;
                        p.P_biased=nan;
                        %p.DeltaS_DG=nan;
                        %P.DeltaS_min=nan;
                        %p.DeltaS_perturb=nan;
                        p.Accuracy=nan;
                    end
                    
                    
                    % p.started=nan;
                    % p.finished=nan;
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
    basename='~/CopyMPI/MaxEnt/Data/PlugInEstimation_1';
    load([basename,'/HQ']);
	numrepeats_save=4;
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
    fprintf('\ndimension: %02d, mean: %02d corr: %02d, std: %02d\n',p.dim, p.fmeanorig, p.startcorr, p.scatterstd)
    if exist([p.filename,'.mat'],'file')
        %try loading it
        A=load([p.filename]);
        results=A.results;
        %keyboard
        if numel(results(1).details)==5 %in this case, means we have old file, so have to do some cleaning up:
            fprintf('\n Doing some cleaning up...')
			for ii=1:numel(results)
				results(ii).details=results(ii).details(1:4);
				results(ii).ZerosInMean=sparse(results(ii).ZerosInMean);
				results(ii).ExitFlag=uint8(results(ii).ExitFlag);
				results(ii).EntropyLog2=[];
				results(ii).Bias_PlugIn=[];
			end
		%keyboard	    		
		save(p.filename,'p','results','fitoptions')    
	elseif numel(results(end).ExitFlag)==numrepeats;
            fprintf('\n File already exists and finished, skipping...')
            continue
		elseif docleanup==false
            fprintf('\n File started, but too lazy to cleanup, moving on\n');
            continue
        else
            kkstart=numel(results(end).ExitFlag)+1;
            fprintf('\n File already exists, starting at repeat %d\n', kkstart)
        end
    else
        kkstart=1;
    end
    
    for kk=kkstart:stoprep
        if mod(kk,50)==1
            tic
            fprintf('Repeat %03d of %d  ', kk, numrepeats);
        end
        
        randn('seed',kk+100*k+1000*i); rand('seed',kk+100*k+1000*i);
        
        featuressampling_all=SampleDiscrete(features,p.P_DG,max(Ns));
        Sigma_q2_true=wcov(features,p.P_q2);
        
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
            
            [lambdalearned,logZlearned, P_q2_learned, meanslearned,output]=FitMaxEntLinear(features,meanssamplingtruncated, fitoptions);
            P_q2_learned=exp(P_q2_learned);
            ExitFlag=output.exitflag;
            Accuracy=max(abs(meanslearned-meanssampling));
            ZerosInMean=sum(meanssampling==0);
            

            EntropyLogE=ent((P_q2_learned));
            EntropyLog2=[]; %EntropyLogE*log2(exp(1));
            
            %now, also have to deal with plug-in estimate: what this means
            %is that we have to calculate Sigma_q2 from the maxent fit to the sampled data,
            %and Sigma_emp_DG by averaging the sampled data, and then to
            %compute their product:
            %keyboard
            Sigma_q2_emp=wcov(features,P_q2_learned);
            Sigma_DG_emp=cov(featuressampling,1);
            %results(k).eig_Sigma_q2_emp=eig(Sigma_q2_emp);
            %results(k).eig_Sigma_DG_emp=eig(Sigma_DG_emp);
warning off
            %Bx=CalcBx(featuressampling,mean(featuressampling),Sigma_q2_emp);
            %Bx_q2_true=CalcBx(featuressampling,p.mean_features_q2,p.Sigma_q2);
            
            if kk<=numrepeats_save
                details.Sigma_q2_emp=Sigma_q2_emp;
                details.Sigma_DG_emp=Sigma_DG_emp;
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
            %results(k).Bias_PlugIn(kk)=mean(Bx);
            results(k).Bias_PlugIn=[]; 
            results(k).Bias_PlugIn_2(kk)=sum(vec(inv(Sigma_q2_emp).*Sigma_DG_emp));
            results(k).Bias_PlugIn_true_q2(kk)=sum(vec(inv(p.Sigma_q2).*Sigma_DG_emp));
            
            results(k).EntropyLogE(kk)=EntropyLogE;
            results(k).EntropyLog2=[];
            results(k).ExitFlag(kk)=int8(ExitFlag);
            results(k).Accuracy(kk)=Accuracy;
            results(k).ZerosInMean(kk)=sparse(ZerosInMean);
            results(k).fitoptions=fitoptions;
           warning on 
       %             keyboard
        end
        %fprintf('\n Done with this round, moving on to next N\n')
        if mod(kk,50)==0
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




