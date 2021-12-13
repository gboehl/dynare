function [resids, grad, state_out, E, M_, out] = match_function(err_0, obs_list,current_obs, opts_simul,...
                                                                            M_, oo_, options_)
% function [resids, grad, stateout, E, M_, out] = match_function(err_0, obs_list,current_obs, opts_simul,...
%                                                                             M_, oo_, options_)
% Outputs:
%  - resids         [double]        [n_exo by 1] vector of residuals
%  - grad           [double]        [n by n_exo] gradient (response of observables to shocks)
%  - state_out      [double]        [ny by 1] value of endogenous variables
%  - E              [double]        response of endogenous variables to shocks
%  - M_             [structure]     Matlab's structure describing the model (M_).
%  - out            [structure]     Occbin's results structure
%
% Inputs
% - err_            [double]        value of shocks 
% - obs_list        [cell]          names of observables
% - current_obs     [double]        [1 by n_obs] current value of observables
% - opts_simul      [structure]     Structure with simulation options
% - M_              [structure]     Matlab's structure describing the model (M_).
% - oo_             [structure]     Matlab's structure containing the results (oo_).
% - options_        [structure]     Matlab's structure describing the current options (options_).

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

opts_simul.SHOCKS = err_0';
options_.occbin.simul=opts_simul;
options_.occbin.simul.full_output=1;
options_.noprint = 1;
[~, out, ss] = occbin.solver(M_,oo_,options_);

nobs = size(obs_list,1);
resids = zeros(nobs,1);

if ~out.error_flag
    state_out= out.piecewise(1,:)' - out.ys;
    
    E = ss.R(:,opts_simul.exo_pos);
    grad = ss.R(opts_simul.varobs_id,opts_simul.exo_pos);
    resids = (out.piecewise(1,opts_simul.varobs_id)-current_obs)'; %-out.endo_ss.(obs_list{this_obs});
else
    grad = NaN(length(opts_simul.varobs_id),length(opts_simul.exo_pos));
    resids = resids+100;
end

end
