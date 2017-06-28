addpath(genpath(pwd))
addpath(genpath('C:\Users\Loki\Desktop\Project\Funs\Misc')) % for minFunc 
%addpath([pwd,'/demo'])
%addpath([pwd,'/temp_files'])
%addpath([pwd,'/maxent'])
%addpath([pwd,'/maxent_MCMC'])
%addpath([pwd,'/dich_gauss'])
%addpath([pwd,'/dich_gauss_bivar_bayes'])
%addpath([pwd,'/flat_models'])
%addpath([pwd,'/util'])
%addpath([pwd,'/third_party'])
%addpath([pwd,'/third_party/minFunc_2012'])
%addpath([pwd,'/third_party/minFunc_2012/minFunc'])
%addpath([pwd,'/third_party/minFunc_2012/minFunc/compiled'])
%addpath([pwd,'/third_party/minFunc_2012/autoDif'])



%compile C_Code  if neceessary
compile_c=0;
if compile_c
    warning('Compiling mex-files, this might take a while');
    cd ./maxent_MCMC/C_Code
    ! make all
    cd ../..
else
    warning('Not compiling mex-files-- set compile_c=1 in file startup.m to compile');
end

