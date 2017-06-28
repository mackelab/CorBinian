%for flat models as well as for "normal" models, want to maximize the bias
%given Delta S-- i.e. to do the opposite of what we usually do (minimize
%Delta S given the bias). To do this, we "simply" invert the function
%MinDeltaS=f(Bias), using matlab-builtin functions, or hand-coded bijection
%search

clear all
addpath([pwd,'/../../FitMaxEnt']);
addpath ~/CopyMPI/Code/minFunc/;
addpath ~/CopyMPI/Code/minFunc/;
addpath ~/CopyMPI/MaxEnt/matlab/MaxEntBias/functions
addpath ~/CopyMPI/BinCorr/matlab/FlatModels
addpath ~/CopyMPI/BinCorr/matlab/FlatDG
addpath ~/CopyMPI/BinCorr/matlab/functions/
addpath ~/CopyMPI/BinCorr/matlab/DGfitting/
basename='~/CopyMPI/MaxEnt/Data/ScatteredMomentsDG_2';
thisdir=pwd;
cd('~/CopyMPI/matlab/functions'), pathchoices=4; startup, cd(thisdir)
minmean=1e-20; %minimum allowed mean, any smaller mean is trunated at this value (otherwise maxent model is not defined, or lives in lower-dimensional space)
fitoptions.TolX=1e-50;
fitoptions.TolFun=1e-50;
fitoptions.display='off';
fitoptions.MaxIter=5000;
fitoptions.restarts=2;
dimo=5;
[features,description,x]=SetupFeaturesMaxEnt(dimo,2);
scalers=linspace(0,1,21);

%first, "normal" models:

gamma=-1;
lambda=.8;

%first, calculate Ising model,
%then, calculate DG model
%then, calculate max bias model at same Delta S
%then, calcualte all minDeltaS models, and compare....


p.gamma=ones(dimo,1)*gamma;
p.Lambda=eye(dimo)*(1-lambda)+ones(dimo)*lambda;
[mu, Rho]=Lambda2Rho(p.gamma, p.Lambda);
p.meanorig=mu(1);
p.corrorig=Rho(2);
p.dim=dimo;

p.P_DG=ProbsDG(x,p.gamma,p.Lambda)';
mean_features_DG=features'*p.P_DG;
p.Entropy_DG=ent(p.P_DG);
[p.lambda_q2,P.logZ_q2, P_q2, p.mean_features_q2,output]=FitMaxEntLinear(features,mean_features_DG', fitoptions);
P_q2=exp(P_q2);
p.Entropy_q2=ent(P_q2);

Sigma_q2=wcov(features,P_q2);
p.Bx=CalcBx(features, p.mean_features_q2,Sigma_q2);
p.Bias_q2=p.Bx'*P_q2;
p.Bias_DG=p.Bx'*p.P_DG;
featuresBiased=[features,p.Bx];
[p.bDash,Hx,Bx,VarB,EBdeltag,SigmaBdeltag]=CalcbDashAtZero(features,p.mean_features_q2, Sigma_q2, P_q2);

[p.Bias_max,p.DeltaS_max, p.lambda_max,p.logZ_max, p.P_max, p.means_max,iterKL,iterBias,output]=FindMaxBiasModel(features, p.Bx, p.mean_features_q2, p.Entropy_q2-p.Entropy_DG, 1e-10, fitoptions)

%keyboard

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

tic
ops.DeltaS=p.Entropy_q2-p.Entropy_DG;
p.model{4}=SampleFromFlatModel('maxBias',p.meanorig,p.corrorig,p.dim,1,1,ops);
[p.model{4}]=CalcBxEtcFlatModel(p.model{4}, p.model{1});
p.model{4}.S=[];
toc

%then lots of minDeltaS models:
for i=1:numel(scalers)
    i
    
    Biases(i)=p.Bias_q2+(p.Bias_DG-p.Bias_q2)*scalers(i);
    ops.Bias=Biases(i);
    p2{i}.model=SampleFromFlatModel('minDeltaS',p.meanorig,p.corrorig,p.dim,1,1,ops);
    [p2{i}.model]=CalcBxEtcFlatModel(p.model{3}, p.model{1});
    

    meansBiased=[p.mean_features_q2,Biases(i)];
    [p2{i}.lambda_biased,p2{i}.logZ_biased, p2{i}.P_biased, p2{i}.means_biased,output]=FitMaxEntLinear(featuresBiased,meansBiased, fitoptions);
    p2{i}.P_biased=exp(p2{i}.P_biased);
    p2{i}.Accuracy=max(abs(p2{i}.means_biased-meansBiased));
    p2{i}.Entropy_biased=ent(p2{i}.P_biased);
    Entropies(i)=p2{i}.Entropy_biased;
    KL(i)=p.Entropy_q2-Entropies(i);
    DeltaSPerturb(i)=(Biases(i)-p.Bias_q2)^2/p.bDash/2;
    KLflat(i)=p2{i}.model.Bias-p.Bias_q2;
end

%%
close all
plot(Biases,KL,'.-');
hold on
plot(Biases,KL,'go-');
plot(Biases,DeltaSPerturb,'-');
plot(p.Bias_DG, p.Entropy_q2-p.Entropy_DG,'kx','linewidth',3);
plot(p.Bias_max, p.DeltaS_max,'ko','linewidth',3);

xlabel('Bias')
ylabel('Delta S')



