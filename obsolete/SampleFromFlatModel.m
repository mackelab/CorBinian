function [model]=SampleFromFlatModel(typer,pspike,rho,n,nsamples,nreps,ops)
warning('Probably better to use CompFlatModels instead')
model.starttime=now;
model.H=nan;
model.H_ana=nan;
model.s=[];
Cpair=rho*pspike*(1-pspike);
pv=repmat(pspike,n,1);
C=repmat(Cpair,n,n)+eye(n)*(pspike*(1-pspike)-Cpair);

r=[0:n];

model.type=typer;
model.pspike=pspike;

model.rho=rho;
model.nreps=nreps;
model.nsamples=nsamples;
model.C=Cpair;
model.exitflag=true;
model.n=n;
model.mucount_target=n*model.pspike;
model.varcount_target=pspike*(1-pspike)*n*(1+(n-1)*rho);

if nargin<=6 || ~isfield(ops,'verbose')
    ops.verbose=1;
end

switch typer
    case 'ising'
       [model.h,model.J,beta,model.Z,model.p_ana,model.exitflag]=FitFlatIsingModel(pspike,Cpair,n);
        [model.HCount_ana,model.H_ana]=EntropyFromSpikeCount(model.p_ana);
        
        if any(~isreal(model.p_ana))
            keyboard
        end
        
    case 'minDeltaS'
        model_q2=model;
        %to get maximally bias model, first need Ising model 
         [model_q2.h,model_q2.J,beta,model_q2.Z,model_q2.p_ana,model_q2.exitflag,junk,output]=FitFlatIsingModel(pspike,Cpair,n);
        [model_q2.mucount_ana,model_q2.varcount_ana]=MeanVar(model_q2.p_ana);
        [model_q2.HCount_ana,model_q2.H_ana]=EntropyFromSpikeCount(model_q2.p_ana);
        %then, need to calculate Bx:
        
        [junk,model_q2]=CalcBxEtcFlatModel(model_q2, model_q2);
        
        [model.h,model.J,model.beta,model.Z,model.p_ana,model.exitflag,junk,output]=FitFlatIsingModel(pspike,Cpair,n,nan,nan,model_q2.Bx_canonical,ops.Bias,model_q2.h,model_q2.J,0);
        [model.HCount_ana,model.H_ana]=EntropyFromSpikeCount(model.p_ana);
        model.model_q2=model_q2;
        model.Bias=model_q2.Bx_canonical'*model.p_ana(:);
        BiasError=abs(model.Bias-ops.Bias);
       % if BiasError>0
         %   keyboard
       % end
        counter=0;
        while counter<=10 && (BiasError>1e-8*(ops.Bias-model_q2.Bias));
            counter;
            counter=counter+1;
            oldmodel=model;
     
        [model.h,model.J,model.beta,model.Z,model.p_ana,model.exitflag,junk,output]=FitFlatIsingModel(pspike,Cpair,n,nan,nan,model_q2.Bx_canonical,ops.Bias,model.h,model.J,model.beta);
        [model.HCount_ana,model.H_ana]=EntropyFromSpikeCount(model.p_ana);
        model.model_q2=model_q2;
        model.Bias=model_q2.Bx_canonical'*model.p_ana(:);
        BiasError=abs(model.Bias-ops.Bias);
           %keyboard 
        end

    case 'maxBias'
        model_q2=model;
        [model_q2.h,model_q2.J,beta,model_q2.Z,model_q2.p_ana,model_q2.exitflag]=FitFlatIsingModel(pspike,Cpair,n);
        [model_q2.mucount_ana,model_q2.varcount_ana]=MeanVar(model_q2.p_ana);
        [model_q2.HCount_ana,model_q2.H_ana]=EntropyFromSpikeCount(model_q2.p_ana);
        [junk,model_q2]=CalcBxEtcFlatModel(model_q2, model_q2);
        model_q2.Cpair=Cpair;
        [model,output]=FindMaxBiasModelFlat(model_q2, ops.DeltaS, 1e-10, ops);
        [model.HCount_ana,model.H_ana]=EntropyFromSpikeCount(model.p_ana);
        model.model_q2=model_q2;
       model.Bias=model_q2.Bx_canonical'*model.p_ana(:);
%       
        
    case 'DG'
       % disp('fitting DG')
        [gamma Lambda] = findLatentGaussian(pv(1:2),C(1:2,1:2));
        gamma=repmat(gamma(1),n,1);
        Lambda=repmat(Lambda(2),n,n)+(Lambda(1)-Lambda(2))*eye(n);
%        keyboard
        [s gamma Lambda Lambdachol] = sampleDichGauss01(gamma,Lambda,1,1,1e-20);
        model.gamma=gamma(1);
        model.lambda=Lambda(2);
        model.p_ana=FlatDGProbs(model.gamma,model.lambda,r,n,5000);
        model.p_ana=model.p_ana/sum(model.p_ana);
        [model.HCount_ana,model.H_ana]=EntropyFromSpikeCount(model.p_ana);
        
   % case 'gibbs'
       % disp('using gibbs sampling')
   %     [model.h,model.J,beta,model.Z,model.p_ana,model.exitflag]=FitFlatIsingModel(pspike,Cpair,n);
   %     model.h_tam=repmat(model.h+model.J/2,n,1);
   %     model.J_tam=triu(repmat(2*model.J,n,n),1);
    case 'niebur'
       % disp('using niebur model')
end

[model.mucount_ana,model.varcount_ana]=MeanVar(model.p_ana);
model.mucount_error=model.mucount_ana-model.mucount_target;
model.varcount_error=model.varcount_ana-model.varcount_target;
model.error=max(abs([model.mucount_error, model.varcount_error]));
switch typer
    case 'maxBias'
        model.DeltaS=(model_q2.H_ana-model.H_ana)/log2(exp(1));
        model.DeltaS_error=abs(ops.DeltaS-model.DeltaS);
        model.DeltaS_desired=ops.DeltaS;
        model.error=max(abs([model.error, model.DeltaS_error]));
        if model.DeltaS_error/ops.DeltaS>1e-3 && ops.DeltaS>1e-5
            disp('DeltaS error too big')
            
            %keyboard
        end
    case 'minDeltaS'
        model.Bias_error=abs(ops.Bias-model.Bias);
        model.error=max(abs([model.error, model.Bias_error]));
        
end

if model.error>.01
    warning('fitting might be inaccurate')
%    keyboard
end



if ~model.exitflag
    model.string='giving up, fitting the model did not work';
    model.H=nan;
    model.mu_acc=nan;
    model.var_acc=nan;
    return
end


%disp('sampling')
model.p=0;
model.mu=0;
model.cov=0;
model.mu_errors=zeros(nreps,1);
model.cov_errors=zeros(nreps,1);
cutter=eye(n)==0;


for k=1:nreps
    switch typer
        case {'ising','maxBias','minDeltaS'}
            [model.s]=SampleFlatIsingSmart(model.p_ana,nsamples);
        case 'DG'
            [model.s gamma junk] = sampleDichGauss01(gamma,Lambda,nsamples,1);
        case 'gibbs'
            %PrintStar(k)
            A=[datestr(clock),'-',num2str(rand,3)];
            A(A=='.' | A==' ' | A==':' | A=='-')='';
            params.tempdir=['/local_jakob/Temp/',A,'/'];
            params.maxtime=60*60;
            model.burninlength=100*numel(model.h_tam);
            model.takeevery=10*numel(model.h_tam);
            %keyboard
            [model.s]=SampleIsing(model.h_tam,model.J_tam,nsamples*model.takeevery+model.burninlength,params)';
            
            model.s=model.s(:,model.burninlength+1:model.takeevery:end);
            
           % keyboard
        case 'niebur'
            model.s= sampleNiebur(pv,Cpair,nsamples,.5)';
    end
    mu=mean(model.s');
    model.mu_errors(k)=sqrt(mean((mu'-pv).^2));
    model.mu=model.mu+mu/nreps;
    c=cov(model.s');
%    model.cov_errors(k)=sqrt(mean((c(cutter)-C(cutter)).^2));
    model.cov=model.cov+c/nreps;
    [p]=CountOnes(model.s');
    [mucounts(k),varcounts(k)]=MeanVar(p);
    model.p=model.p+p;
end
model.p=model.p/sum(model.p);
[model.mucount_samp,model.varcount_samp]=MeanVar(model.p);
model.mu_count_mean=mean(mucounts);
model.mu_count_std=std(mucounts);
model.var_count_mean=mean(varcounts);
model.var_count_std=std(varcounts);




[model.Hcount,model.H]=EntropyFromSpikeCount(model.p);
s=model.type;
s=[s,repmat(' ',1,6-numel(s))];
if ops.verbose
    if isfield(model,'H_ana') H=model.H_ana; else H=model.H; end
%fprintf(['\n %s  entropy  %3f mu_acc %3f var_acc %3f'],s, H,mean(model.mu_errors),mean(model.cov_errors))
end
model.time=(now-model.starttime)*24*60*60;
%keyboard

