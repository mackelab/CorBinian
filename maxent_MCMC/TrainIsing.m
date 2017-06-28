function [results,params]=TrainIsing(mu_true, cov_true, params)
%function [results,params]=TrainIsing(mu_true, cov_true, params)
%
%fit Ising model
%
%outputs:
%results.h  are the biases (local fields)
%results.J are the couplings
%results.track is lots of info about the fitting
%params has all the parameters used for fitting

%inputs:
%mu_true, cov_true the desired means and covariances
%params is a struct with fitting options, leave empty for default
%
%metrics: l1 distance on means, l1 distance on covariances, l2 distance on
%correlation coefficients
%
%
%co Jakob H Macke, Jakob.Macke@gmail.com 2010
%using some methods originally by T Broderick and others.
% see http://arxiv.org/pdf/0712.2437.pdf
% and http://www.ncbi.nlm.nih.gov/pubmed/22539825


%params
%(if any field is left empty, then the default entry is used)
%.tempdir: Directory in which to save intermediate results
%.innerls=10, Number of inner-loop iterations in optimization
%.outerls=1000, Number of outer loops
%.MC=10000, size of each MCMC sample
%.finalMC=10000, size of final MCMC sample
%.maxtime=60*60*24, maximum runtime of algorithm in seconds
%.exit_conditions: Three numbers indicating when 'convergence' is achieved.
%Default is using function ApproxFinishLines to set them heuristically
%.remove_files=1, whether to remove temporary files at the end
%.collect_all_params, whether to also store the parameters after each outer
%loop
%.verbose=1, whether to print out stuff or not
%.exit_condition_scaler=1: A number by which all numbers in exit_conditions
%are scaled
%.re_train=0: Whther to restart the algorithm if it does not convere
%.init='init': Initialize parameters by indep model

if nargin==3
    p=params;
else
    p=struct;
end


%make a temporary directory in which to dump all the files:
results.starttime=now;

p.mu_true=mu_true(:)';
p.cov_true=cov_true;

BaseDir=pwd;
%BaseLoc=[BaseDir 'Retina/Shapes/matlab/Decoders/Ising/'];


%disp('First, have to set the two environment variables: Directory where the IsingCode is')
BaseLoc=[base_dir '/maxent_MCMC/'];

%disp('Directory where temporary results should be saved')
TempLoc=[base_dir '/temp_files/'];
%keyboard




rand('seed',randn);
A=[datestr(clock),'-',num2str(rand,3)];
A(A=='.' | A==' ' | A==':' | A=='-')='';

%set a couple of default parameters, in fact, quite many....
%if ~isfield(p,'tempdir') p.tempdir=[BaseLoc,'temp_files/',A,'/']; end
if ~isfield(p,'tempdir') p.tempdir=[TempLoc,A,'/']; end

mkdir(p.tempdir);
%keyboard
if ~isfield(p,'innerls') p.innerls=10; end  % number of of inner-loop iterations
if ~isfield(p,'outerls') p.outerls=1000; end %numer of outer-loop iterations
if ~isfield(p,'MC') p.MC=10000; end %size of each MC sample
if ~isfield(p,'finalMC') p.finalMC=p.MC; end   %size of last MC-sampe
if ~isfield(p,'maxtime') p.maxtime=60*60*24; end %maximal runtime of algorithm (well, it will its current loop before exiting)
if ~isfield(p,'exit_conditions') p.exit_conditions=ApproxFinishLines(p.MC,mu_true,cov_true); end
if ~isfield(p,'remove_files') p.remove_files=1; end
if ~isfield(p,'collect_all_params') p.collect_all_params=1; end
if ~isfield(p,'flip_it') p.flip_it=1; end
if ~isfield(p,'verbose') p.verbose=1; end
if ~isfield(p,'exit_condition_scaler'), p.exit_condition_scaler=1; end
if ~isfield(p,'re_train'), p.re_train=0; end
if ~isfield(p,'init') p.init='indep'; end

p.exit_conditions=p.exit_conditions.*p.exit_condition_scaler;

%set names of temporary files
p.TCODE=[BaseLoc 'C_Code/'];
p.dataFile=[p.tempdir,'Cors.txt'];
p.paramsFile=[p.tempdir,'learnedParams'];
p.trackFile=[p.tempdir,'TrackFile.txt'];
p.logFile=[p.tempdir,'LogFile.txt'];
p.outputFile=[p.tempdir,'output.txt'];
p.MonteFile=[p.tempdir,'Monte.txt'];


ncells = length(p.mu_true);

p.ncells=ncells;
if p.maxtime<0 p.maxtime=abs(p.maxtime)*p.ncells; end

if p.flip_it
    p.flipper=double(p.mu_true>.5);
    %p.flipper=0*p.mu_true;
    %p.flipper(1:5)=1;
else
    p.flipper=zeros(size(p.mu_true));
end
flipmatrix=diag(-2*p.flipper+1);
p.mu_train=p.mu_true*flipmatrix+p.flipper;
p.cov_train=flipmatrix*p.cov_true*flipmatrix;
p.secmom_train=p.cov_train+p.mu_train'*p.mu_train;

WriteDataFile(p.dataFile,p.mu_train,p.secmom_train);
%keyboard

%if no intial conditions are specified, take independent model

%keyboard
switch p.init
    case 'indep'
        p.initCommand='indep';
    case 'writetofile'
        p.init='write initial params to file';
        if ~isfield(p,'initFile') p.initFile=[p.tempdir,'initParams.txt']; end
        p.initCommand=p.initFile;
        %    keyboard
        J2=flipmatrix*p.Jinit*flipmatrix;
        h2=flipmatrix*(p.hinit+.5*(J2+J2')*p.flipper');
        %    keyboard
        WriteDataFile(p.initFile,h2,J2/2);
    case 'findinfile'
        p.initCommand=p.initFile;
    otherwise 
        error
end




save([p.tempdir,'params'],'p','results')

%make a big command string out of all the options
p.inputs=['-s_' num2str(p.MC) ' -i_' num2str(p.innerls) ' -n_' num2str(p.ncells) ...
    ' -d_' num2str(p.finalMC)    ' -r_' num2str(p.outerls) ' -J_' p.initCommand ...
    ' -t_' num2str(p.maxtime) ' -u_' num2str(p.exit_conditions(1)) ' -v_' num2str(p.exit_conditions(2)) ' -w_' num2str(p.exit_conditions(3))];

%make an even bigger comand string
p.command=[p.TCODE 'learn_stop.exe ' p.inputs ' -D_' p.dataFile ' -L_' p.trackFile ' -P_' p.paramsFile ' -M_' p.MonteFile ' &> ' p.logFile];

if p.verbose
    disp(['Starting Ising model fitting procedure ', datestr(now)])
    disp(['using ' num2str(p.ncells) ' dimensions.'])
    disp(['Be patient, this can take some time!'])
    disp(['Progress can be monitored via log-files in ' p.tempdir])
    disp(['Finish lines are ', num2str(p.exit_conditions)]);
end


%keyboard
[p.status p.output]=system(p.command);
%keyboard

%read final parameter estimates
[results.h,results.J]=ReadResultDataFile([p.paramsFile,'.txt'],p.flipper);

%read all the diagnoistics into results.track
%keyboard
results.track = ReadTrackFile(p.trackFile);

if (p.collect_all_params || all(isnan(results.h))) && results.track.success

    %read parameters after each outer loop
    for k=1:max(results.track.outer)
        [results.hprelim(:,k),results.Jprelim(:,:,k)]=ReadResultDataFile([p.paramsFile,num2str(k),'.txt'],p.flipper);
    end

    %read out empirical means after each outer loop
    for k=1:max(results.track.outer)
        [muprelim,covprelim]=ReadResultDataFile([p.trackFile,'Moments',num2str(k),'.txt']);
        covprelim=covprelim/2;
        covprelim=covprelim+covprelim'+diag(muprelim);
        covprelim=covprelim-muprelim(:)*muprelim(:)';
        muprelim=muprelim'*flipmatrix+p.flipper;
        covprelim=flipmatrix*covprelim*flipmatrix;
        results.cov_prelim(:,:,k)=covprelim;
        results.mu_prelim(:,k)=muprelim;
    end
end

results.log=fscanf(fopen(p.logFile,'r'),'%s');

%if results.h is nan, the result-datfile could not be found. In this case,
%read of the last of the preliminary results, but issue a warning.
if all(isnan(results.h))
    results.status='Final results could not be read, use preliminary results instead';
    warning(results.status)
    try
        results.h=results.hprelim(:,end);
        results.J=results.Jprelim(:,:,end);
        if p.collect_all_params==0
            results.hprelim=[];
            results.Jprelim=[];
        end
    catch
        results.status='preliminary results could not be read';
        warning(results.status);
    end
end


%read out last MC sample, if it exists;
try
    Monte = load(p.MonteFile);
    ncells = Monte(1);
    Monte = Monte(2:end);
    Monte = reshape(Monte, length(Monte)/ncells, ncells);
    Monte=Monte*flipmatrix+repmat(p.flipper,size(Monte,1),1);
    results.mu = mean(Monte, 1);
    results.cov = cov(Monte);
    results.MC=logical(Monte);
    %results.MC_Pchange=mixtrace(results.MC);
catch
    warning('Error loading last MC sample')
    results.MC='last MC sample could not be found';
end


if p.ncells<=15
    [results.Z,results.states,results.mu_hJ,results.cov_hJ,results.allPs]=CompZ(results.h,results.J);
else
    results.Z=1;
    %[results.Z,junk,results.mu_hJ,results.cov_hJ,junk]=CompZ(results.h,results.J,results.MC);

end

results.mu_true=p.mu_true;
results.cov_true=p.cov_true;
results.endtime=now;
results.traintime=(results.endtime-results.starttime)*24*60;
params=p;


%delete temporary file
if p.remove_files
    try
        delete([p.tempdir,'*']);
        rmdir(p.tempdir)
    catch
        warning('Problems removing temporary files')
        %   keyboard
    end
end




if p.verbose
    try
        disp(['Optimization finished because: ',results.track.exit]);
        disp(['Optimization took ', num2str(results.traintime) ' minutes'])
        disp(['Errors on means, covariances, corr-coeffs:'])
        fprintf('%1.1E %1.1E %1.1E \n',results.track.metrics(end,:));
        disp(['Finish lines were:'])
        fprintf('%1.1E %1.1E %1.1E \n',p.exit_conditions);
        disp([''])
    catch
        warning('track file could not be read')
    end
end

%keyboard
if ~results.track.converged && params.re_train>0
    params.re_train=p.re_train-1;
    params.MC=params.MC*2;
    cov_scaled=p.regparameter*cov_true;
    cov_scaled=cov_scaled-diag(diag(cov_scaled))+diag(diag(cov_true));
    disp(['Not converged, starting again with twice as many samples, and ',num2str(p.regparameter),' of correlations'])
    disp(['After this, I will try ', num2str(params.re_train),' more times.'])
    [results,params]=TrainIsing(mu_true, cov_scaled, params);

end


function A=mixtrace(samples)
A=abs(diff(samples));
A=mean(A,2);


function A=Cov2Corr(A)
%takes in a covariance matrix,  returns the corresponding correlation
%matrix
%std(A)=diag(A)
stdA=repmat(sqrt(diag(A)),1,size(A,1));
A=A./stdA./(stdA');

