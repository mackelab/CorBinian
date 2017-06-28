function [model]=CompFlatModels(typer,pspike,rho,n,nsamples,ops)
model.starttime=now;
model.H=nan;
model.s=[];

Cpair=rho*pspike*(1-pspike);
pv=repmat(pspike,n,1);
C=repmat(Cpair,2,2)+eye(2)*(pspike*(1-pspike)-Cpair);

r=[0:n];

model.type=typer;
model.pspike=pspike;
model.rho=rho;
model.nsamples=nsamples;
model.C=Cpair;
model.CC=C;
model.exitflag=true;
model.n=n;


if nargin<=5
    ops.verbose=1;
    ops.acc=1e-10;
    ops.steps=5000;
    ops.logbase=2;
    ops.verbose=0;
end
if ~isfield(ops,'acc') ops.acc=1e-100; end
if ~isfield(ops,'steps') ops.steps=5000; end
if ~isfield(ops,'verbose') ops.verbose=0; end


switch typer
    case 'indep'
        model.p=binopdf(r,n,pspike);
    case 'ising'
        % disp('fitting Ising model')
        [model.h,model.J,beta,model.Z,model.p,model.exitflag,model.lognchoosek]=FitFlatIsingModel(pspike,Cpair,n,ops.acc);
    case 'ising3'
        [model.h,model.J,model.gamma,beta,model.Z,model.p,model.exitflag,model.lognchoosek]=FitFlatIsingModelThird(pspike,Cpair,ops.E3count,n,ops.acc);    
    case 'isinghJ'
        model.h=pspike;
        model.J=rho;
        [model.p,model.Z]=IsiSpikeCountDistrib(model.h,model.J,n);
    case 'minDeltaS'
        model_q2=model;
        %to get maximally bias model, first need Ising model 
         [model_q2.h,model_q2.J,beta,model_q2.Z,model_q2.p,model_q2.exitflag,junk,output]=FitFlatIsingModel(pspike,Cpair,n);

        [junk,model_q2]=CalcBxEtcFlatModel(model_q2, model_q2);
        
        [model.h,model.J,model.beta,model.Z,model.p_ana,model.exitflag,junk,output]=FitFlatIsingModel(pspike,Cpair,n,nan,nan,model_q2.Bx_canonical,ops.Bias,model_q2.h,model_q2.J,0);
        model.model_q2=model_q2;
        model.Bias=model_q2.Bx_canonical'*model.p(:);
        BiasError=abs(model.Bias-ops.Bias);
     
        counter=0;
        while counter<=10 && (BiasError>1e-8*(ops.Bias-model_q2.Bias));
            counter;
            counter=counter+1;
            oldmodel=model;
     
        [model.h,model.J,model.beta,model.Z,model.p,model.exitflag,junk,output]=FitFlatIsingModel(pspike,Cpair,n,nan,nan,model_q2.Bx_canonical,ops.Bias,model.h,model.J,model.beta);
        model.model_q2=model_q2;
        model.Bias=model_q2.Bx_canonical'*model.p(:);
        BiasError=abs(model.Bias-ops.Bias);
           %keyboard 
        end

    case 'maxBias'
        model_q2=model;
        [model_q2.h,model_q2.J,beta,model_q2.Z,model_q2.p,model_q2.exitflag]=FitFlatIsingModel(pspike,Cpair,n);

        [junk,model_q2]=CalcBxEtcFlatModel(model_q2, model_q2);
        model_q2.Cpair=Cpair;
        model.model_q2=model_q2;
       model.Bias=model_q2.Bx_canonical'*model.p(:);
       
    case 'DG'
        % disp('fitting DG')
        if isfield(ops,'gamma');
            model.gamma=ops.gamma;
            model.lambda=ops.lambda;
        else
%            keyboard
        C=repmat(Cpair,2,2)+eye(2)*(pspike*(1-pspike)-Cpair);
        pv=repmat(pspike,2,1);
        [gamma Lambda] = findLatentGaussian(pv(1:2),C(1:2,1:2),max(1e-5,ops.acc));
        model.gamma=gamma(1);
        model.lambda=Lambda(2);
        end
      %  keyboard
        model.p=FlatDGProbs(model.gamma,model.lambda,r,n,ops.steps);
end

if isfield(ops,'beta')
    %use temperature beta instead of 1
    beta=ops.beta;
    lognchoosek=(gammaln(n+1)-gammaln(r+1)-gammaln(n-r+1));
    logp=(1-ops.beta)*lognchoosek+beta*log(model.p);
    model.logZ=logsumexp(logp');
    model.Z=exp(model.logZ);
    model.p=exp(logp-model.logZ);
end

model.p=model.p/sum(model.p);
[model.mucount,model.varcount]=MeanVar(model.p);
 pspike=model.mucount/n;
 model.pspike=pspike;
 model.varspike=pspike*(1-pspike);
 model.C=(model.varcount-n*model.varspike)/n/(n-1);
 model.rho=model.C/model.varspike;  
[model.HCount,model.H]=EntropyFromSpikeCount(model.p,ops.logbase);
[model.HC]=HeatCapacityFromSpikeCount(model.p,ops.logbase);


if sum(isnan(model.p))>0 || isnan(model.H) || model.H<.00001
  %  keyboard
end


if ~model.exitflag
    model.string='giving up, fitting the model did not work';
    model.H=nan;
    model.mu_acc=nan;
    model.var_acc=nan;
    return
end



[model.s]=SampleFlatIsingSmart(model.p,nsamples);



if ops.verbose
    s=model.type;
    s=[s,repmat(' ',1,6-numel(s))];
    fprintf(['\n %s  entropy  %3f'],s, model.H)
end

model.time=(now-model.starttime)*24*60*60;

