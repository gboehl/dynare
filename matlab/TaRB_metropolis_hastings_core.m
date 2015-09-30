function myoutput = TaRB_metropolis_hastings_core(myinputs,fblck,nblck,whoiam, ThisMatlab)
% function myoutput = TaRB_metropolis_hastings_core(myinputs,fblck,nblck,whoiam, ThisMatlab)
% Contains the most computationally intensive portion of code in
% random_walk_metropolis_hastings (the 'for xxx = fblck:nblck' loop) using the TaRB algorithm. 
% The branches in  that 'for'
% cycle are completely independent to be suitable for parallel execution.
%
% INPUTS
%   o myimput            [struc]     The mandatory variables for local/remote
%                                    parallel computing obtained from random_walk_metropolis_hastings.m
%                                    function.
%   o fblck and nblck    [integer]   The Metropolis-Hastings chains.
%   o whoiam             [integer]   In concurrent programming a modality to refer to the different threads running in parallel is needed.
%                                    The integer whoaim is the integer that
%                                    allows us to distinguish between them. Then it is the index number of this CPU among all CPUs in the
%                                    cluster.
%   o ThisMatlab         [integer]   Allows us to distinguish between the
%                                    'main' Matlab, the slave Matlab worker, local Matlab, remote Matlab,
%                                     ... Then it is the index number of this slave machine in the cluster.
% OUTPUTS
%   o myoutput  [struc]
%               If executed without parallel, this is the original output of 'for b =
%               fblck:nblck'. Otherwise, it's a portion of it computed on a specific core or
%               remote machine. In this case:
%                               record;
%                               irun;
%                               NewFile;
%                               OutputFileName
%
% ALGORITHM
%   Portion of Tailored Randomized Block Metropolis-Hastings proposed in
%   Chib/Ramamurthy (2010): Tailored randomized block MCMC methods with
%   application to DSGE models, Journal of Econometrics 155, pp. 19-38
% 
%   This implementation differs from the originally proposed one in the
%   treatment of non-positive definite Hessians. Here we
%       - use the Jordan decomposition
%
% SPECIAL REQUIREMENTS.
%   None.
% 
% PARALLEL CONTEXT
% See the comments in the random_walk_metropolis_hastings.m funtion.


% Copyright (C) 2006-2015 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if nargin<4,
    whoiam=0;
end

% reshape 'myinputs' for local computation.
% In order to avoid confusion in the name space, the instruction struct2local(myinputs) is replaced by:

TargetFun=myinputs.TargetFun;
ProposalFun=myinputs.ProposalFun;
xparam1=myinputs.xparam1;
mh_bounds=myinputs.mh_bounds;
last_draw=myinputs.ix2;
last_posterior=myinputs.ilogpo2;
fline=myinputs.fline;
npar=myinputs.npar;
nruns=myinputs.nruns;
NewFile=myinputs.NewFile;
MAX_nruns=myinputs.MAX_nruns;
d=myinputs.d;
InitSizeArray=myinputs.InitSizeArray;
record=myinputs.record;
dataset_ = myinputs.dataset_;
dataset_info = myinputs.dataset_info;
bayestopt_ = myinputs.bayestopt_;
estim_params_ = myinputs.estim_params_;
options_ = myinputs.options_;
M_ = myinputs.M_;
oo_ = myinputs.oo_;

% Necessary only for remote computing!
if whoiam
    % initialize persistent variables in priordens()
    priordens(xparam1,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7, bayestopt_.p3,bayestopt_.p4,1);
end

MetropolisFolder = CheckPath('metropolis',M_.dname);
ModelName = M_.fname;
BaseName = [MetropolisFolder filesep ModelName];

options_.lik_algo = 1;
OpenOldFile = ones(nblck,1);

%
% Now I run the (nblck-fblck+1) Metropolis-Hastings chains
%

block_iter=0;

for curr_chain = fblck:nblck,
    block_iter=block_iter+1;
    try
        % This will not work if the master uses a random number generator not
        % available in the slave (different Matlab version or
        % Matlab/Octave cluster). Therefore the trap.
        %
        % Set the random number generator type (the seed is useless but needed by the function)
        set_dynare_seed(options_.DynareRandomStreams.algo, options_.DynareRandomStreams.seed);
        % Set the state of the RNG
        set_dynare_random_generator_state(record.InitialSeeds(curr_chain).Unifor, record.InitialSeeds(curr_chain).Normal);
    catch
        % If the state set by master is incompatible with the slave, we only reseed 
        set_dynare_seed(options_.DynareRandomStreams.seed+curr_chain);
    end
    if (options_.load_mh_file~=0) && (fline(curr_chain)>1) && OpenOldFile(curr_chain) %load previous draws and likelihood
        load([BaseName '_mh' int2str(NewFile(curr_chain)) '_blck' int2str(curr_chain) '.mat'])
        x2 = [x2;zeros(InitSizeArray(curr_chain)-fline(curr_chain)+1,npar)];
        logpo2 = [logpo2;zeros(InitSizeArray(curr_chain)-fline(curr_chain)+1,1)];
        OpenOldFile(curr_chain) = 0;
    else
        x2 = zeros(InitSizeArray(curr_chain),npar);
        logpo2 = zeros(InitSizeArray(curr_chain),1);
    end
    %Prepare waiting bars
    if whoiam
        prc0=(curr_chain-fblck)/(nblck-fblck+1)*(isoctave || options_.console_mode);
        hh = dyn_waitbar({prc0,whoiam,options_.parallel(ThisMatlab)},['MH (' int2str(curr_chain) '/' int2str(options_.mh_nblck) ')...']);
    else
        hh = dyn_waitbar(0,['Metropolis-Hastings (' int2str(curr_chain) '/' int2str(options_.mh_nblck) ')...']);
        set(hh,'Name','Metropolis-Hastings');
    end
    accepted_draws_this_chain = 0;
    accepted_draws_this_file = 0;
    blocked_draws_counter_this_chain=0;
    blocked_draws_counter_this_chain_this_file=0;
    draw_index_current_file = fline(curr_chain); %get location of first draw in current block
    draw_iter = 1;
    while draw_iter <= nruns(curr_chain)
        
        %% randomize indices for blocking in this iteration
        indices=randperm(npar)';
        blocks=[1; (1+cumsum((rand(length(indices)-1,1)>(1-options_.TaRB.new_block_probability))))];
        nblocks=blocks(end,1); %get number of blocks this iteration
        current_draw=last_draw(curr_chain,:)'; %get starting point for current draw for updating
        
        for block_iter=1:nblocks
            blocked_draws_counter_this_chain=blocked_draws_counter_this_chain+1;
            blocked_draws_counter_this_chain_this_file=blocked_draws_counter_this_chain_this_file+1;
            nxopt=length(indices(blocks==block_iter,1)); %get size of current block
            par_start_current_block=current_draw(indices(blocks==block_iter,1));
            [xopt_current_block, fval, exitflag, hess_mat_optimizer, options_, Scale] = dynare_minimize_objective(@TaRB_optimizer_wrapper,par_start_current_block,options_.TaRB.mode_compute,options_,[mh_bounds.lb(indices(blocks==block_iter,1),1) mh_bounds.ub(indices(blocks==block_iter,1),1)],bayestopt_.name,bayestopt_,[],...
                current_draw,indices(blocks==block_iter,1),TargetFun,...% inputs for wrapper
                dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_); %inputs for objective           
            %% covariance for proposal density
            hessian_mat = reshape(hessian('TaRB_optimizer_wrapper',xopt_current_block, ...
                    options_.gstep,...
                    current_draw,indices(blocks==block_iter,1),TargetFun,...% inputs for wrapper
                dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_),nxopt,nxopt);
            
            if any(any(isnan(hessian_mat))) || any(any(isinf(hessian_mat)))
                inverse_hessian_mat=eye(nxopt)*1e-4; %use diagonal
            else
                inverse_hessian_mat=inv(hessian_mat); %get inverse Hessian
                if any(any((isnan(inverse_hessian_mat)))) || any(any((isinf(inverse_hessian_mat))))
                    inverse_hessian_mat=eye(nxopt)*1e-4; %use diagonal
                end               
            end
            [proposal_covariance_Cholesky_decomposition_upper,negeigenvalues]=chol(inverse_hessian_mat);            
            %if not positive definite, use generalized Cholesky if
            %Eskow/Schnabel
            if negeigenvalues~=0 
                proposal_covariance_Cholesky_decomposition_upper=chol_SE(inverse_hessian_mat,0);
            end
            proposal_covariance_Cholesky_decomposition_upper=proposal_covariance_Cholesky_decomposition_upper*diag(bayestopt_.jscale(indices(blocks==block_iter,1),:));
            %get proposal draw
            if strcmpi(ProposalFun,'rand_multivariate_normal')
                n = nxopt;
            elseif strcmpi(ProposalFun,'rand_multivariate_student')
                n = options_.student_degrees_of_freedom;
            end
    
            proposed_par = feval(ProposalFun, xopt_current_block', proposal_covariance_Cholesky_decomposition_upper, n);
            % chech whether draw is valid and compute posterior
            if all( proposed_par(:) > mh_bounds.lb(indices(blocks==block_iter,1),:) ) && all( proposed_par(:) < mh_bounds.ub(indices(blocks==block_iter,1),:) )
                try
                    logpost = - feval('TaRB_optimizer_wrapper', proposed_par(:),...
                        current_draw,indices(blocks==block_iter,1),TargetFun,...% inputs for wrapper
                        dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_);
                catch
                    logpost = -inf;
                end
            else
                logpost = -inf;
            end
            %get ratio of proposal densities, required because proposal depends
            %on current mode via Hessian and is thus not symmetric anymore
            if strcmpi(ProposalFun,'rand_multivariate_normal')
                proposal_density_proposed_move_forward=multivariate_normal_pdf(proposed_par,xopt_current_block',proposal_covariance_Cholesky_decomposition_upper,n);
                proposal_density_proposed_move_backward=multivariate_normal_pdf(par_start_current_block',xopt_current_block',proposal_covariance_Cholesky_decomposition_upper,n);
            elseif strcmpi(ProposalFun,'rand_multivariate_student')
                proposal_density_proposed_move_forward=multivariate_student_pdf(proposed_par,xopt_current_block',proposal_covariance_Cholesky_decomposition_upper,n);
                proposal_density_proposed_move_backward=multivariate_student_pdf(par_start_current_block',xopt_current_block',proposal_covariance_Cholesky_decomposition_upper,n);
            end
            accprob=logpost-last_posterior(curr_chain)+ log(proposal_density_proposed_move_backward)-log(proposal_density_proposed_move_forward); %Formula (6), Chib/Ramamurthy

            if (logpost > -inf) && (log(rand) < accprob)
                current_draw(indices(blocks==block_iter,1))=proposed_par;
                last_posterior(curr_chain)=logpost;
                accepted_draws_this_chain =accepted_draws_this_chain +1;
                accepted_draws_this_file = accepted_draws_this_file + 1;
            else %no updating 
                %do nothing, keep old value
            end
        end
        %save draws and update stored last values
        x2(draw_index_current_file,:) = current_draw;
        last_draw(curr_chain,:) = current_draw;
        %save posterior after full run through all blocks        
        logpo2(draw_index_current_file) = last_posterior(curr_chain);

        prtfrc = draw_iter/nruns(curr_chain);
        if (mod(draw_iter, 3)==0 && ~whoiam) || (mod(draw_iter,50)==0 && whoiam)
            dyn_waitbar(prtfrc,hh,[ 'MH (' int2str(curr_chain) '/' int2str(options_.mh_nblck) ') ' sprintf('Current acceptance ratio %4.3f', accepted_draws_this_chain/blocked_draws_counter_this_chain)]);
        end
        if (draw_index_current_file == InitSizeArray(curr_chain)) || (draw_iter == nruns(curr_chain)) % Now I save the simulations, either because the current file is full or the chain is done
            save([BaseName '_mh' int2str(NewFile(curr_chain)) '_blck' int2str(curr_chain) '.mat'],'x2','logpo2');
            fidlog = fopen([MetropolisFolder '/metropolis.log'],'a');
            fprintf(fidlog,['\n']);
            fprintf(fidlog,['%% Mh' int2str(NewFile(curr_chain)) 'Blck' int2str(curr_chain) ' (' datestr(now,0) ')\n']);
            fprintf(fidlog,' \n');
            fprintf(fidlog,['  Number of simulations.: ' int2str(length(logpo2)) '\n']);
            fprintf(fidlog,['  Acceptance ratio......: ' num2str(accepted_draws_this_file /blocked_draws_counter_this_chain_this_file) '\n']);
            fprintf(fidlog,['  Posterior mean........:\n']);
            for i=1:length(x2(1,:))
                fprintf(fidlog,['    params:' int2str(i) ': ' num2str(mean(x2(:,i))) '\n']);
            end
            fprintf(fidlog,['    log2po:' num2str(mean(logpo2)) '\n']);
            fprintf(fidlog,['  Minimum value.........:\n']);
            for i=1:length(x2(1,:))
                fprintf(fidlog,['    params:' int2str(i) ': ' num2str(min(x2(:,i))) '\n']);
            end
            fprintf(fidlog,['    log2po:' num2str(min(logpo2)) '\n']);
            fprintf(fidlog,['  Maximum value.........:\n']);
            for i=1:length(x2(1,:))
                fprintf(fidlog,['    params:' int2str(i) ': ' num2str(max(x2(:,i))) '\n']);
            end
            fprintf(fidlog,['    log2po:' num2str(max(logpo2)) '\n']);
            fprintf(fidlog,' \n');
            fclose(fidlog);
            %reset counters;
            accepted_draws_this_file = 0;
            blocked_draws_counter_this_chain_this_file=0;
            if draw_iter == nruns(curr_chain) % I record the last draw...
                record.LastParameters(curr_chain,:) = x2(end,:);
                record.LastLogPost(curr_chain) = logpo2(end);
            end
            % size of next file in chain curr_chain
            InitSizeArray(curr_chain) = min(nruns(curr_chain)-draw_iter,MAX_nruns);
            % initialization of next file if necessary
            if InitSizeArray(curr_chain)
                x2 = zeros(InitSizeArray(curr_chain),npar);
                logpo2 = zeros(InitSizeArray(curr_chain),1);
                NewFile(curr_chain) = NewFile(curr_chain) + 1;
                draw_index_current_file = 0;
            end
        end
        draw_iter=draw_iter+1;
        draw_index_current_file = draw_index_current_file + 1;
    end% End of the simulations for one mh-block.
    record.AcceptanceRatio(curr_chain) = accepted_draws_this_chain/blocked_draws_counter_this_chain;
    dyn_waitbar_close(hh);
    [record.LastSeeds(curr_chain).Unifor, record.LastSeeds(curr_chain).Normal] = get_dynare_random_generator_state();
    OutputFileName(block_iter,:) = {[MetropolisFolder,filesep], [ModelName '_mh*_blck' int2str(curr_chain) '.mat']};
end% End of the loop over the mh-blocks.


myoutput.record = record;
myoutput.irun = draw_index_current_file;
myoutput.NewFile = NewFile;
myoutput.OutputFileName = OutputFileName;