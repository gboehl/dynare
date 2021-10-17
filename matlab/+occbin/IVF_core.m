function [filtered_errs, resids, Emat, stateval, error_code] = IVF_core(M_,oo_,options_,err_index,filtered_errs_init,my_obs_list,obs,init_val)
% function [filtered_errs, resids, Emat, stateval] = IVF_core(M_,oo_,options_,err_index,filtered_errs_init,my_obs_list,obs,init_val)
% Computes thre 
%
% Outputs:
%  - filtered_errs          [T by N_obs]        filtered shocks
%  - resids                 [T by N_obs]        residuals
%  - Emat                   [N by N_obs by T]   response matrix of endogenous variables to shocks at each point in time
%  - stateval               [T by N]            vector of endogenous variables
%  - error_code             [4 by 1]            error code 
%
% Inputs
% - M_                      [structure]     Matlab's structure describing the model (M_).
% - oo_                     [structure]     Matlab's structure containing the results (oo_).
% - options_                [structure]     Matlab's structure describing the current options (options_).
% - err_index               [double]        index of shocks with strictly positive variance in M_.exo_names
% - filtered_errs_init      [T by N_obs]    initial values for the shocks
% - my_obs_list             [cell]          names of observables
% - obs                     [T by N_obs]    observed data        
% - init_val                [N by 1]        initial value of endogenous variables

% Original authors: Pablo Cuba-Borda, Luca Guerrieri, Matteo Iacoviello, and Molin Zhong
% Original file downloaded from:
% http://www.lguerrieri.com/jae-replication.zip
% Adapted for Dynare by Dynare Team.
%
% This code is in the public domain and may be used freely.
% However the authors would appreciate acknowledgement of the source by
% citation of any of the following papers:
%
% Pablo Cuba-Borda, Luca Guerrieri, Matteo Iacoviello, and Molin Zhong (2019): "Likelihood evaluation of models 
% with occasionally binding constraints", Journal of Applied Econometrics,
% 34(7), 1073-1085


%-------------------------------------
% Filter shocks
%-------------------------------------

[sample_length, n_obs]= size(obs);
nerrs = size(err_index,1);
if nargin<8
    init_val = zeros(M_.endo_nbr,1);
end

resids = zeros(sample_length,nerrs);
stateval = zeros(sample_length,M_.endo_nbr);
Emat = zeros(M_.endo_nbr,nerrs,sample_length);
error_code = zeros(4,1);
%solver options (set locally)
options_.solve_algo = options_.occbin.solver.solve_algo;
options_.solve_tolf = options_.occbin.solver.solve_tolf;
options_.solve_tolx = options_.occbin.solver.solve_tolx;
options_.options.steady.maxit = options_.occbin.solver.maxit;
options_.jacobian_flag=1;

opts_simul = options_.occbin.simul;
opts_simul.curb_retrench = options_.occbin.likelihood.curb_retrench;
opts_simul.maxit = options_.occbin.likelihood.maxit;
opts_simul.waitbar = false;
opts_simul.periods = options_.occbin.likelihood.periods;
opts_simul.check_ahead_periods = options_.occbin.likelihood.check_ahead_periods;
opts_simul.periodic_solution = options_.occbin.likelihood.periodic_solution;
opts_simul.restrict_state_space = options_.occbin.likelihood.restrict_state_space;
opts_simul.piecewise_only = 1;

filtered_errs=zeros(sample_length,n_obs);

if options_.occbin.likelihood.waitbar
    hh = dyn_waitbar(0,'IVF_core: Filtering the shocks');
    set(hh,'Name','IVF_core: Filtering the shocks.');
end
    
for this_period=1:sample_length
    if options_.occbin.likelihood.waitbar
        dyn_waitbar(this_period/sample_length, hh, sprintf('Period %u of %u', this_period,sample_length));
    end
    current_obs = obs(this_period,:);
    init_val_old = init_val;
    
    inan = ~isnan(current_obs);
    current_obs = current_obs(inan);
    obs_list = my_obs_list(inan);
    opts_simul.varobs_id=options_.varobs_id(inan)';
    opts_simul.exo_pos=err_index(inan); %err_index is predefined mapping from observables to shocks
    opts_simul.endo_init = init_val_old;
    opts_simul.SHOCKS = filtered_errs_init(this_period,inan);
%         [ err_vals_out, exitflag ] = csolve(@(err_vals) occbin.match_function(...
%             err_vals, obs_list,current_obs, opts_simul, M_,oo_,options_),...
%             err0',@(err_vals) occbin.match_function(...
%             err_vals, obs_list,current_obs, opts_simul, M_,oo_,options_),options_.solve_tolf,options_.occbin.solver.maxit);
    [err_vals_out, exitflag] = dynare_solve(@occbin.match_function, filtered_errs_init(this_period,inan)', options_, obs_list,current_obs, opts_simul, M_,oo_,options_);
    
    if exitflag
        filtered_errs=NaN;
        error_code(1) = 304;
        error_code(4) = 1000;
        if options_.occbin.likelihood.waitbar; dyn_waitbar_close(hh); end
        return    
    end
    filtered_errs(this_period,inan)=err_vals_out';
    
    opts_simul.SHOCKS = err_vals_out;
    
    [ resids(this_period,inan), ~, stateval(this_period,:), Emat(:,inan,this_period), M_] = occbin.match_function(...
        err_vals_out,obs_list,current_obs,opts_simul, M_,oo_,options_);
    init_val = stateval(this_period,:); %update
    if max(abs(err_vals_out))>1e8
        error_code(1) = 306;
        error_code(4) = max(abs(err_vals_out))/1000;
        filtered_errs=NaN;
        if options_.occbin.likelihood.waitbar; dyn_waitbar_close(hh); end
        return
    end
    if max(abs(resids(this_period,:)))>0.001
        disp_verbose('IVF_core: I am stopping because match_function could not find the shocks that',options_.verbosity)
        disp_verbose('IVF_core: solve for the model''s observed variables',options_.verbosity)
        filtered_errs=NaN;
        error_code(1) = 303;
        error_code(4) = max(abs(resids(this_period,:)))*100;
        if options_.occbin.likelihood.waitbar; dyn_waitbar_close(hh); end
        return
    end
end
if options_.occbin.likelihood.waitbar
    dyn_waitbar_close(hh); 
end

end