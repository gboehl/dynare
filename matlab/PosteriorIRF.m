function PosteriorIRF(type)
% Builds posterior IRFs after the MH algorithm.
%
% INPUTS
%   o type       [char]     string specifying the joint density of the
%                           deep parameters ('prior','posterior').
%
% OUTPUTS
%   None                    (oo_ and plots).
%
% SPECIAL REQUIREMENTS
%   None

% PARALLEL CONTEXT
% This funtion has been parallelized in two different points. Then we have two core
% functions associated with it(the _core1 and _core2).
% See also the comments posterior_sampler.m funtion.

% Copyright © 2006-2018 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.


global options_ estim_params_ oo_ M_ bayestopt_ dataset_ dataset_info

% Set the number of periods
if isempty(options_.irf) || ~options_.irf
    options_.irf = 40;
end
% Set varlist if necessary
varlist = options_.varlist;
if isempty(varlist)
    varlist = options_.varobs;
end
options_.varlist = varlist;
nvar = length(varlist);
IndxVariables = [];
for i=1:nvar
    idx = strmatch(varlist{i}, M_.endo_names, 'exact');
    if isempty(idx)
        disp(['PosteriorIRF :: ' varlist{i} 'is not a declared endogenous variable!'])
    else
        IndxVariables = [IndxVariables, idx];
    end
end

% Get index of shocks for requested IRFs
irf_shocks_indx = getIrfShocksIndx();

% Set various parameters & Check or create directories
nvx  = estim_params_.nvx;
nvn  = estim_params_.nvn;
ncx  = estim_params_.ncx;
ncn  = estim_params_.ncn;
np   = estim_params_.np ;
npar = nvx+nvn+ncx+ncn+np;
offset = npar-np; clear('nvx','nvn','ncx','ncn','np');

nvobs = dataset_.vobs;
gend = dataset_.nobs;
MaxNumberOfPlotPerFigure = 9;
nn = sqrt(MaxNumberOfPlotPerFigure);
MAX_nirfs_dsge = ceil(options_.MaxNumberOfBytes/(options_.irf*nvar*M_.exo_nbr)/8);
MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);
if options_.dsge_var
    MAX_nirfs_dsgevar = ceil(options_.MaxNumberOfBytes/(options_.irf*nvobs*M_.exo_nbr)/8);
else
    MAX_nirfs_dsgevar = 0;
end

DirectoryName = CheckPath('Output',M_.dname);
if strcmpi(type,'posterior')
    MhDirectoryName = CheckPath('metropolis',M_.dname);
elseif strcmpi(type,'gsa')
    if options_.opt_gsa.pprior
        MhDirectoryName = CheckPath(['GSA' filesep 'prior'],M_.dname);
    else
        MhDirectoryName = CheckPath(['GSA' filesep 'mc'],M_.dname);
    end
else
    MhDirectoryName = CheckPath('prior',M_.dname);
end

%delete old stale files before creating new ones
delete_stale_file([MhDirectoryName filesep M_.fname '_IRF_DSGEs*.mat']);
delete_stale_file([MhDirectoryName filesep M_.fname '_IRF_BVARDSGEs*.mat']);
delete_stale_file([MhDirectoryName filesep M_.fname '_irf_dsge*.mat']);
delete_stale_file([MhDirectoryName filesep M_.fname '_irf_bvardsge*.mat']);
delete_stale_file([MhDirectoryName filesep M_.fname '_param_irf*.mat']);


if strcmpi(type,'posterior')
    B = options_.sub_draws;
    options_.B = B;
    if round((1-options_.mh_conf_sig)*B)<2
        fprintf('\nPosteriorIRF:: options_.mh_conf_sig times options_.sub_draws is too small to generate HPDIs. I am omitting them.\n')
    end
elseif strcmpi(type,'gsa')
    RootDirectoryName = CheckPath('gsa',M_.dname);
    if options_.opt_gsa.pprior
        load([ RootDirectoryName filesep  M_.fname '_prior.mat'],'lpmat0','lpmat','istable')
    else
        load([ RootDirectoryName filesep  M_.fname '_mc.mat'],'lpmat0','lpmat','istable')
    end
    x=[lpmat0(istable,:) lpmat(istable,:)];
    clear lpmat istable
    B=size(x,1); options_.B = B;
else% type = 'prior'
    B = options_.prior_draws;
    options_.B = B;
end

irun = 0;
IRUN = 0;
irun2 = 0;
NumberOfIRFfiles_dsge = 1;
NumberOfIRFfiles_dsgevar = 1;
ifil2 = 1;
% Create arrays
if B <= MAX_nruns
    stock_param = zeros(B, npar);
else
    stock_param = zeros(MAX_nruns, npar);
end
if B >= MAX_nirfs_dsge
    stock_irf_dsge = zeros(options_.irf,nvar,M_.exo_nbr,MAX_nirfs_dsge);
else
    stock_irf_dsge = zeros(options_.irf,nvar,M_.exo_nbr,B);
end
if MAX_nirfs_dsgevar
    if B >= MAX_nirfs_dsgevar
        stock_irf_bvardsge = zeros(options_.irf,nvobs,M_.exo_nbr,MAX_nirfs_dsgevar);
    else
        stock_irf_bvardsge = zeros(options_.irf,nvobs,M_.exo_nbr,B);
    end
    NumberOfLags = options_.dsge_varlag;
    NumberOfLagsTimesNvobs = NumberOfLags*nvobs;
    if options_.noconstant
        NumberOfParametersPerEquation = NumberOfLagsTimesNvobs;
    else
        NumberOfParametersPerEquation = NumberOfLagsTimesNvobs+1;
    end
    Companion_matrix = diag(ones(nvobs*(NumberOfLags-1),1),-nvobs);
end

% First block of code executed in parallel. The function devoted to do it is  PosteriorIRF_core1.m
% function.

b = 0;

localVars=[];

% Save the local variables.

localVars.IRUN = IRUN;
localVars.irun = irun;
localVars.irun2=irun2;
localVars.npar = npar;

localVars.type=type;
if strcmpi(type,'posterior')
    while b<B
        b = b + 1;
        x(b,:) = GetOneDraw(type,M_,estim_params_,oo_,options_,bayestopt_);
    end
end

if ~strcmpi(type,'prior')
    localVars.x=x;
end

b=0;
if options_.dsge_var
    localVars.gend = gend;
    localVars.nvobs = nvobs;
    localVars.NumberOfParametersPerEquation = NumberOfParametersPerEquation;
    localVars.NumberOfLags = options_.dsge_varlag;
    localVars.NumberOfLagsTimesNvobs = NumberOfLags*nvobs;
    localVars.Companion_matrix = diag(ones(nvobs*(NumberOfLags-1),1),-nvobs);
end
localVars.nvar=nvar;
localVars.IndxVariables=IndxVariables;
localVars.MAX_nirfs_dsgevar=MAX_nirfs_dsgevar;
localVars.MAX_nirfs_dsge=MAX_nirfs_dsge;
localVars.MAX_nruns=MAX_nruns;
localVars.NumberOfIRFfiles_dsge=NumberOfIRFfiles_dsge;
localVars.NumberOfIRFfiles_dsgevar=NumberOfIRFfiles_dsgevar;
localVars.ifil2=ifil2;
localVars.MhDirectoryName=MhDirectoryName;

% Like sequential execution!
if isnumeric(options_.parallel)
    [fout] = PosteriorIRF_core1(localVars,1,B,0);
    nosaddle = fout.nosaddle;
else
    % Parallel execution!
    [nCPU, totCPU, nBlockPerCPU] = distributeJobs(options_.parallel, 1, B);
    for j=1:totCPU-1
        nfiles = ceil(nBlockPerCPU(j)/MAX_nirfs_dsge);
        NumberOfIRFfiles_dsge(j+1) =NumberOfIRFfiles_dsge(j)+nfiles;
        if MAX_nirfs_dsgevar
            nfiles = ceil(nBlockPerCPU(j)/MAX_nirfs_dsgevar);
        else
            nfiles=0;
        end
        NumberOfIRFfiles_dsgevar(j+1) =NumberOfIRFfiles_dsgevar(j)+nfiles;
        nfiles = ceil(nBlockPerCPU(j)/MAX_nruns);
        ifil2(j+1) =ifil2(j)+nfiles;
    end
    localVars.NumberOfIRFfiles_dsge=NumberOfIRFfiles_dsge;
    localVars.NumberOfIRFfiles_dsgevar=NumberOfIRFfiles_dsgevar;
    localVars.ifil2=ifil2;

    globalVars = struct('M_',M_, ...
                        'options_', options_, ...
                        'bayestopt_', bayestopt_, ...
                        'estim_params_', estim_params_, ...
                        'oo_', oo_, ...
                        'dataset_',dataset_, ...
                        'dataset_info',dataset_info);

    % which files have to be copied to run remotely
    NamFileInput(1,:) = {'',[M_.fname '.static.m']};
    NamFileInput(2,:) = {'',[M_.fname '.dynamic.m']};
    if M_.set_auxiliary_variables
        NamFileInput(3,:) = {'',[M_.fname '.set_auxiliary_variables.m']};
    end
    if options_.steadystate_flag
        if options_.steadystate_flag == 1
            NamFileInput(length(NamFileInput)+1,:)={'',[M_.fname '_steadystate.m']};
        else
            NamFileInput(length(NamFileInput)+1,:)={'',[M_.fname '.steadystate.m']};
        end
    end
    [fout] = masterParallel(options_.parallel, 1, B,NamFileInput,'PosteriorIRF_core1', localVars, globalVars, options_.parallel_info);
    nosaddle=0;
    for j=1:length(fout)
        nosaddle = nosaddle + fout(j).nosaddle;
    end

end

% END first parallel section!

if nosaddle
    disp(['PosteriorIRF :: Percentage of discarded posterior draws = ' num2str(nosaddle/(B+nosaddle))])
end

ReshapeMatFiles('irf_dsge',type)
if MAX_nirfs_dsgevar
    ReshapeMatFiles('irf_bvardsge')
end

if strcmpi(type,'gsa')
    return
end

IRF_DSGEs = dir([MhDirectoryName filesep M_.fname '_IRF_DSGEs*.mat']);
NumberOfIRFfiles_dsge = length(IRF_DSGEs);

IRF_BVARDSGEs = dir([MhDirectoryName filesep M_.fname '_IRF_BVARDSGEs*.mat']);
NumberOfIRFfiles_dsgevar = length(IRF_BVARDSGEs);

MeanIRF = zeros(options_.irf,nvar,M_.exo_nbr);
MedianIRF = zeros(options_.irf,nvar,M_.exo_nbr);
VarIRF = zeros(options_.irf,nvar,M_.exo_nbr);
DistribIRF = zeros(options_.irf,9,nvar,M_.exo_nbr);
HPDIRF = zeros(options_.irf,2,nvar,M_.exo_nbr);

if options_.TeX
    varlist_TeX = cell(nvar, 1);
    for i=1:nvar
        varlist_TeX(i) = {M_.endo_names_tex{IndxVariables(i)}};
    end
end

fprintf('Estimation::mcmc: Posterior (dsge) IRFs...\n');
tit(M_.exo_names_orig_ord) = M_.exo_names;
kdx = 0;

for file = 1:NumberOfIRFfiles_dsge
    load([MhDirectoryName filesep M_.fname '_IRF_DSGEs' int2str(file) '.mat']);
    for i = irf_shocks_indx
        for j = 1:nvar
            for k = 1:size(STOCK_IRF_DSGE,1)
                kk = k+kdx;
                [MeanIRF(kk,j,i),MedianIRF(kk,j,i),VarIRF(kk,j,i),HPDIRF(kk,:,j,i),...
                 DistribIRF(kk,:,j,i)] = posterior_moments(squeeze(STOCK_IRF_DSGE(k,j,i,:)),0,options_.mh_conf_sig);
            end
        end
    end
    kdx = kdx + size(STOCK_IRF_DSGE,1);

end

clear STOCK_IRF_DSGE;

for i = irf_shocks_indx
    for j = 1:nvar
        name = sprintf('%s_%s', M_.endo_names{IndxVariables(j)}, tit{i});
        oo_.PosteriorIRF.dsge.Mean.(name) = MeanIRF(:,j,i);
        oo_.PosteriorIRF.dsge.Median.(name) = MedianIRF(:,j,i);
        oo_.PosteriorIRF.dsge.Var.(name) = VarIRF(:,j,i);
        oo_.PosteriorIRF.dsge.deciles.(name) = DistribIRF(:,:,j,i);
        oo_.PosteriorIRF.dsge.HPDinf.(name) = HPDIRF(:,1,j,i);
        oo_.PosteriorIRF.dsge.HPDsup.(name) = HPDIRF(:,2,j,i);
    end
end


if MAX_nirfs_dsgevar
    MeanIRFdsgevar = zeros(options_.irf,nvar,M_.exo_nbr);
    MedianIRFdsgevar = zeros(options_.irf,nvar,M_.exo_nbr);
    VarIRFdsgevar = zeros(options_.irf,nvar,M_.exo_nbr);
    DistribIRFdsgevar = zeros(options_.irf,9,nvar,M_.exo_nbr);
    HPDIRFdsgevar = zeros(options_.irf,2,nvar,M_.exo_nbr);
    fprintf('Estimation::mcmc: Posterior (bvar-dsge) IRFs...\n');
    tit(M_.exo_names_orig_ord) = M_.exo_names;
    kdx = 0;
    for file = 1:NumberOfIRFfiles_dsgevar
        load([MhDirectoryName filesep M_.fname '_IRF_BVARDSGEs' int2str(file) '.mat']);
        for i = irf_shocks_indx
            for j = 1:nvar
                for k = 1:size(STOCK_IRF_BVARDSGE,1)
                    kk = k+kdx;
                    [MeanIRFdsgevar(kk,j,i),MedianIRFdsgevar(kk,j,i),VarIRFdsgevar(kk,j,i),...
                     HPDIRFdsgevar(kk,:,j,i),DistribIRFdsgevar(kk,:,j,i)] = ...
                        posterior_moments(squeeze(STOCK_IRF_BVARDSGE(k,j,i,:)),0,options_.mh_conf_sig);
                end
            end
        end
        kdx = kdx + size(STOCK_IRF_BVARDSGE,1);
    end
    clear STOCK_IRF_BVARDSGE;
    for i = irf_shocks_indx
        for j = 1:nvar
            name = sprintf('%s_%s', M_.endo_names{IndxVariables(j)}, tit{i});
            oo_.PosteriorIRF.bvardsge.Mean.(name) = MeanIRFdsgevar(:,j,i);
            oo_.PosteriorIRF.bvardsge.Median.(name) = MedianIRFdsgevar(:,j,i);
            oo_.PosteriorIRF.bvardsge.Var.(name) = VarIRFdsgevar(:,j,i);
            oo_.PosteriorIRF.bvardsge.deciles.(name) = DistribIRFdsgevar(:,:,j,i);
            oo_.PosteriorIRF.bvardsge.HPDinf.(name) = HPDIRFdsgevar(:,1,j,i);
            oo_.PosteriorIRF.bvardsge.HPDsup.(name) = HPDIRFdsgevar(:,2,j,i);
        end
    end
end
%
%      Finally I build the plots.
%


% Second block of code executed in parallel, with the exception of file
% .tex generation always run in sequentially. This portion of code is execute in parallel by
% PosteriorIRF_core2.m function.

if ~options_.nograph && ~options_.no_graph.posterior
    % Save the local variables.
    localVars=[];

    Check=options_.TeX;
    if (Check)
        localVars.varlist_TeX=varlist_TeX;
    end


    localVars.nvar=nvar;
    localVars.MeanIRF=MeanIRF;
    localVars.tit=tit;
    localVars.nn=nn;
    localVars.MAX_nirfs_dsgevar=MAX_nirfs_dsgevar;
    localVars.HPDIRF=HPDIRF;
    localVars.varlist=varlist;
    localVars.MaxNumberOfPlotPerFigure=MaxNumberOfPlotPerFigure;
    if options_.dsge_var
        localVars.HPDIRFdsgevar=HPDIRFdsgevar;
        localVars.MeanIRFdsgevar = MeanIRFdsgevar;
    end

    % The files .TeX are genereted in sequential way always!

    subplotnum = 0;
    titTeX(M_.exo_names_orig_ord) = M_.exo_names_tex;
    if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        fidTeX = fopen([DirectoryName filesep M_.fname '_BayesianIRF.tex'],'w');
        fprintf(fidTeX,'%% TeX eps-loader file generated by PosteriorIRF.m (Dynare).\n');
        fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
        fprintf(fidTeX,' \n');
        titTeX(M_.exo_names_orig_ord) = M_.exo_names_tex;

        for ii=irf_shocks_indx
            figunumber = 0;

            for jj=1:nvar
                if max(abs(MeanIRF(:,jj,ii))) >= options_.impulse_responses.plot_threshold
                    subplotnum = subplotnum+1;

                    if subplotnum == 1
                        fprintf(fidTeX,'\\begin{figure}[H]\n');
                    end
                end
                if subplotnum == MaxNumberOfPlotPerFigure || (jj == nvar  && subplotnum> 0)
                    figunumber = figunumber+1;

                    fprintf(fidTeX,'\\centering \n');
                    fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s/%s_Bayesian_IRF_%s_%d}\n',options_.figures.textwidth*min(subplotnum/nn,1),DirectoryName,M_.fname,tit{ii},figunumber);
                    if options_.relative_irf
                        fprintf(fidTeX,['\\caption{Bayesian relative IRF.}']);
                    else
                        fprintf(fidTeX,'\\caption{Bayesian IRF: Orthogonalized shock to $%s$.}\n',titTeX{ii});
                    end
                    fprintf(fidTeX,'\\label{Fig:BayesianIRF:%s:%d}\n', tit{ii},figunumber);
                    fprintf(fidTeX,'\\end{figure}\n');
                    fprintf(fidTeX,' \n');
                    subplotnum = 0;
                end
            end
        end
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
    end

    % The others file format are generated in parallel by PosteriorIRF_core2!

    if ~isoctave
        if isnumeric(options_.parallel)  || (M_.exo_nbr*ceil(length(varlist)/MaxNumberOfPlotPerFigure))<8
            [fout] = PosteriorIRF_core2(localVars,1,M_.exo_nbr,0);
        else
            isRemoteOctave = 0;
            for indPC=1:length(options_.parallel)
                isRemoteOctave = isRemoteOctave + (findstr(options_.parallel(indPC).MatlabOctavePath, 'octave'));
            end
            if isRemoteOctave
                [fout] = PosteriorIRF_core2(localVars,1,M_.exo_nbr,0);
            else
                globalVars = struct('M_',M_, ...
                                    'options_', options_);

                [fout] = masterParallel(options_.parallel, 1, M_.exo_nbr,NamFileInput,'PosteriorIRF_core2', localVars, globalVars, options_.parallel_info);
            end
        end
    else
        [fout] = PosteriorIRF_core2(localVars,1,M_.exo_nbr,0);
    end
    % END parallel code!

end

fprintf('Estimation::mcmc: Posterior IRFs, done!\n');
