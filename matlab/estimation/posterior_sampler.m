function posterior_sampler(TargetFun,ProposalFun,xparam1,sampler_options,mh_bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_,dispString)
% function posterior_sampler(TargetFun,ProposalFun,xparam1,sampler_options,mh_bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_,dispString)
% Random Walk Metropolis-Hastings algorithm.
%
% INPUTS
%   o TargetFun         [char]      string specifying the name of the objective
%                                   function (posterior kernel).
%   o ProposalFun       [char]      string specifying the name of the proposal
%                                   density
%   o xparam1           [double]    (p*1) vector of parameters to be estimated (initial values).
%   o sampler_options               structure
%   o mh_bounds         [double]    (p*2) matrix defining lower and upper bounds for the parameters.
%   o dataset_          [structure] data structure (likelihood-based) or data moments (method of moments)
%   o dataset_info      [structure] dataset info structure (likelihood-based) or info on weighting matrix (method of moments)
%   o options_          [structure] options structure
%   o M_                [structure] model structure
%   o estim_params_     [structure] estimated parameters structure
%   o bayestopt_        [structure] prior specification structure
%   o oo_               [structure] output structure
%   o dispString        [string]    string prependening the messages printed to the command window
%
% SPECIAL REQUIREMENTS
%   None.
%
% PARALLEL CONTEXT
% The most computationally intensive part of this function may be executed
% in parallel. The code suitable to be executed in
% parallel on multi core or cluster machine (in general a 'for' cycle)
% has been removed from this function and been placed in the posterior_sampler_core.m funtion.
%
% The DYNARE parallel packages comprise a i) set of pairs of Matlab functions that can be executed in
% parallel and called name_function.m and name_function_core.m and ii) a second set of functions used
% to manage the parallel computations.
%
% This function was the first function to be parallelized. Later, other
% functions have been parallelized using the same methodology.
% Then the comments write here can be used for all the other pairs of
% parallel functions and also for management functions.

% Copyright © 2006-2023 Dynare Team
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

vv = sampler_options.invhess;
% Initialization of the sampler
[ ix2, ilogpo2, ModelName, MetropolisFolder, fblck, fline, npar, nblck, nruns, NewFile, MAX_nruns, d, bayestopt_] = ...
    posterior_sampler_initialization(TargetFun, xparam1, vv, mh_bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_, dispString);

InitSizeArray = min([repmat(MAX_nruns,nblck,1) fline+nruns-1],[],2);

% Load last mh history file
record=load_last_mh_history_file(MetropolisFolder, ModelName);

% Only for test parallel results!!!

% To check the equivalence between parallel and serial computation!
% First run in serial mode, and then comment the follow line.
%   save('recordSerial.mat','-struct', 'record');

% For parallel runs after serial runs with the abobe line active.
%   TempRecord=load('recordSerial.mat');
%   record.Seeds=TempRecord.Seeds;



% Snapshot of the current state of computing. It necessary for the parallel
% execution (i.e. to execute in a corretct way a portion of code remotely or
% on many cores). The mandatory variables for local/remote parallel
% computing are stored in the localVars struct.

localVars =   struct('TargetFun', TargetFun, ...
                     'ProposalFun', ProposalFun, ...
                     'xparam1', xparam1, ...
                     'vv', vv, ...
                     'sampler_options', sampler_options, ...
                     'mh_bounds', mh_bounds, ...
                     'ix2', ix2, ...
                     'ilogpo2', ilogpo2, ...
                     'ModelName', ModelName, ...
                     'fline', fline, ...
                     'npar', npar, ...
                     'nruns', nruns, ...
                     'NewFile', NewFile, ...
                     'MAX_nruns', MAX_nruns, ...
                     'd', d, ...
                     'InitSizeArray',InitSizeArray, ...
                     'record', record, ...
                     'dataset_', dataset_, ...
                     'dataset_info', dataset_info, ...
                     'options_', options_, ...
                     'M_',M_, ...
                     'bayestopt_', bayestopt_, ...
                     'estim_params_', estim_params_, ...
                     'dr', oo_.dr,...
                     'endo_steady_state', oo_.steady_state,...
                     'exo_steady_state', oo_.exo_steady_state,...
                     'exo_det_steady_state', oo_.exo_det_steady_state,...
                     'varargin',[]);

if strcmp(sampler_options.posterior_sampling_method,'tailored_random_block_metropolis_hastings')
    localVars.options_.silent_optimizer=1; %locally set optimizer to silent mode
    if ~isempty(sampler_options.optim_opt)
        localVars.options_.optim_opt=sampler_options.optim_opt; %locally set options for optimizer
    end
end

% User doesn't want to use parallel computing, or wants to compute a
% single chain compute sequentially.

if isnumeric(options_.parallel) || (~isempty(fblck) && (nblck-fblck)==0)
    fout = posterior_sampler_core(localVars, fblck, nblck, 0);
    record = fout.record;
    % Parallel in Local or remote machine.
else
    % Global variables for parallel routines.
    globalVars = struct();
    % which files have to be copied to run remotely
    NamFileInput(1,:) = {'',[ModelName '.static.m']};
    NamFileInput(2,:) = {'',[ModelName '.dynamic.m']};
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
    if (options_.load_mh_file~=0)  && any(fline>1)
        NamFileInput(length(NamFileInput)+1,:)={[M_.dname '/metropolis/'],[ModelName '_mh' int2str(NewFile(1)) '_blck*.mat']};
    end
    % from where to get back results
    %     NamFileOutput(1,:) = {[M_.dname,'/metropolis/'],'*.*'};
    if options_.mh_recover && isempty(fblck)
        % here we just need to retrieve the output of the completed remote jobs
        fblck=1;
        options_.parallel_info.parallel_recover = 1;
    end
    [fout, nBlockPerCPU, totCPU] = masterParallel(options_.parallel, fblck, nblck,NamFileInput,'posterior_sampler_core', localVars, globalVars, options_.parallel_info);
    for j=1:totCPU
        offset = sum(nBlockPerCPU(1:j-1))+fblck-1;
        record.LastLogPost(offset+1:sum(nBlockPerCPU(1:j)))=fout(j).record.LastLogPost(offset+1:sum(nBlockPerCPU(1:j)));
        record.LastParameters(offset+1:sum(nBlockPerCPU(1:j)),:)=fout(j).record.LastParameters(offset+1:sum(nBlockPerCPU(1:j)),:);
        record.AcceptanceRatio(offset+1:sum(nBlockPerCPU(1:j)))=fout(j).record.AcceptanceRatio(offset+1:sum(nBlockPerCPU(1:j)));
        record.FunctionEvalPerIteration(offset+1:sum(nBlockPerCPU(1:j)))=fout(j).record.FunctionEvalPerIteration(offset+1:sum(nBlockPerCPU(1:j)));
        record.LastSeeds(offset+1:sum(nBlockPerCPU(1:j)))=fout(j).record.LastSeeds(offset+1:sum(nBlockPerCPU(1:j)));
        if j==1
            if isfield(fout(j).record,'ProposalCovariance') && isfield(fout(j).record,'ProposalScaleVec')
                record.ProposalCovariance=fout(j).record.ProposalCovariance;
                record.ProposalScaleVec=fout(j).record.ProposalScaleVec;
            end
        end
    end
    options_.parallel_info.parallel_recover = 0;
end

irun = fout(1).irun;
NewFile = fout(1).NewFile;

record.MCMCConcludedSuccessfully = 1; %set indicator for successful run

update_last_mh_history_file(MetropolisFolder, ModelName, record);

% Provide diagnostic output
skipline()
fprintf('%s: Number of mh files: %u per block.\n',    dispString, NewFile(1));
fprintf('%s: Total number of generated files: %u.\n', dispString, NewFile(1)*nblck);
fprintf('%s: Total number of iterations: %u.\n',      dispString, (NewFile(1)-1)*MAX_nruns+irun-1);
fprintf('%s: Current acceptance ratio per chain:\n', dispString);
for i=1:nblck
    if i<10
        disp(['                                                       Chain  ' num2str(i) ': ' num2str(100*record.AcceptanceRatio(i)) '%'])
    else
        disp(['                                                       Chain ' num2str(i) ': ' num2str(100*record.AcceptanceRatio(i)) '%'])
    end
end
if max(record.FunctionEvalPerIteration)>1
    fprintf('%s: Current function evaluations per iteration:\n', dispString);
    for i=1:nblck
        if i<10
            disp(['                                                       Chain  ' num2str(i) ': ' num2str(record.FunctionEvalPerIteration(i))])
        else
            disp(['                                                       Chain ' num2str(i) ': ' num2str(record.FunctionEvalPerIteration(i))])
        end
    end
end
