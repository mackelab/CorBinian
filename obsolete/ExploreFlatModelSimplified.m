
%% %% %%

%simulate from Ising model and DG for increasing population size, and fit
%bias analytically and numerically in both cases, to see how the bias
%scales with population size:

if 1
    clear all
    addpath([pwd,'/../../FitMaxEnt']);
    addpath ~/CopyMPI/Code/minFunc/;
    addpath ~/CopyMPI/MaxEnt/matlab/MaxEntBias/functions
    addpath ~/CopyMPI/BinCorr/matlab/FlatModels
    addpath ~/CopyMPI/BinCorr/matlab/FlatDG
    addpath ~/CopyMPI/BinCorr/matlab/functions/
    addpath ~/CopyMPI/BinCorr/matlab/DGfitting/
    basename='~/CopyMPI/MaxEnt/Data/ExloreFlatModelSimplified';
    
    mkdir(basename)
    callfile=mfilename;
end


dims=14; %use Ising model with 2, 5, 10, 15 dimensions
%dims=5;

means=.2 %try out means 0.01, .1, .5
rhos=.5; %take pairwise correlation 0, .1, .5

%options for parameter learning:
fitoptions.TolX=1e-50;
fitoptions.TolFun=1e-50;
fitoptions.display='off';
fitoptions.MaxIter=3000;


%first, set up all the parameters and save them:
%    savenameparams=[basename,'Params'];
randn('seed',0); rand('seed',0);


p=struct;
p.dim=dims;
p.meanorig=means;
p.corrorig=rhos;
[features,description,x]=SetupFeaturesMaxEnt(p.dim,2);

%first, do stuff for infinite range Ising model:
%warning off
[model]=SampleFromFlatModel('ising',p.meanorig,p.corrorig,p.dim,1,1);
p.PCount_q2=model.p_ana;
counts=sum(x,2);
p.P_q2=p.PCount_q2(counts+1)';
for i=1:numel(p.P_q2);
    p.P_q2(i)=p.P_q2(i)/nchoosek(p.dim,counts(i));
end

%now, get the mean from the _q2 model, and fit _q2
%model to it (cumbersome, but want to have same
%numerical error everywhere, and too lazy to think
%about converting parameters...)
p.means_q2=p.P_q2'*features;
p.Entropy_q2=ent(p.P_q2);


[p.lambda_q2,p.logZ_q2, junk, p.means_q2,output]=FitMaxEntLinear(features,p.means_q2, fitoptions);

model_q2=model;
[model]=SampleFromFlatModel('DG',p.meanorig,p.corrorig,p.dim,1,1);
%warning on
p.PCount_DG=model.p_ana;
%counts=sum(x,2);
p.P_DG=p.PCount_DG(counts+1)';
for i=1:numel(p.P_DG);
    p.P_DG(i)=p.P_DG(i)/nchoosek(p.dim,counts(i));
end
p.means_DG=p.P_DG'*features;
[p.lambda_DG,p.logZ_DG, junk, p.means_DG,output]=FitMaxEntLinear(features,p.means_DG, fitoptions);
p.Entropy_DG=ent(p.P_DG);

p.Sigma_q2=wcov(features,p.P_q2);
p.Bx=CalcBx(features, p.means_q2,p.Sigma_q2);



p.Bias_q2=p.Bx'*p.P_q2;
p.Bias_DG=p.Bx'*p.P_DG;
ops.Bias=p.Bias_DG;
[model_maxbias]=SampleFromFlatModel('maxbias',p.meanorig,p.corrorig,p.dim,1,1,ops);
p.PCount_maxbias=model_maxbias.p_ana;
p.P_maxbias=p.PCount_maxbias(counts+1)';
for i=1:numel(p.P_maxbias);
    p.P_maxbias(i)=p.P_maxbias(i)/nchoosek(p.dim,counts(i));
end


p.Bias_maxbias=p.Bx'*p.P_maxbias;


%set up features for max-ent fitting with constraints on bias:
featuresBiased=[features,p.Bx];
[p.bDash,Hx,Bx,VarB,EBdeltag,SigmaBdeltag]=CalcbDashAtZero(features,p.means_q2, p.Sigma_q2, p.P_q2);

%now, calculate distribution which has minimal
%higher-order correlations, but same bias as DG:
meansBiased=[p.means_q2,p.Bias_DG];
try
    [p.lambda_biased,p.logZ_biased, p.P_biased, p.means_biased,output]=FitMaxEntLinear(featuresBiased,meansBiased, fitoptions);
    p.P_biased=exp(p.P_biased);
    p.Accuracy=max(abs(p.means_biased-meansBiased));
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

p.Entropy_biased=ent(p.P_biased);
p.DeltaS_DG=p.Entropy_q2-p.Entropy_DG;
p.DeltaS_min=p.Entropy_q2-p.Entropy_biased;
p.DeltaS_perturb=(p.Bias_q2-p.Bias_DG).^2/p.bDash;
p.Entropy_perturb=p.Entropy_q2-p.DeltaS_perturb;

[model,model_q2]=CalcBxEtcFlatModel(model, model_q2)

%%
%at the moment: 
%find maximal 

%plots that I want, as function of parameters mu and rho:
%
%additional bias of DG
%additional bias of maximally biased model
%perturbative prediction of maximally biased model
%bdash as function of parameters







