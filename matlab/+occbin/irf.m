function irfs = irf(M_,oo_,options_)
% irfs = irf(M_,oo_,options_)
% Calls a minimizer
%
% INPUTS
% - M_                  [structure]     Matlab's structure describing the model
% - oo_                 [structure]     Matlab's structure containing the results
% - options_            [structure]     Matlab's structure describing the current options
%
% OUTPUTS
% - irfs                [structure]     IRF results
%
% SPECIAL REQUIREMENTS
%   none.
%
%
% Copyright © 2022-2023 Dynare Team
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

shocknames = options_.occbin.irf.exo_names;
shocksigns = options_.occbin.irf.shocksigns; %'pos','neg'
shocksize = options_.occbin.irf.shocksize;
t_0 = options_.occbin.irf.t0;

%% set simulation options based on IRF options
options_.occbin.simul.init_regime = options_.occbin.irf.init_regime;
options_.occbin.simul.check_ahead_periods = options_.occbin.irf.check_ahead_periods;
options_.occbin.simul.maxit = options_.occbin.irf.maxit;
options_.occbin.simul.periods = options_.irf;

%% Run initial conditions + other shocks
if t_0 == 0
    shocks_base = zeros(options_.occbin.simul.periods+1,M_.exo_nbr);
    options_.occbin.simul.endo_init = [];
else
    if ~isfield(oo_.occbin,'smoother')
        error('occbin.irfs: smoother must be run before requesting GIRFs based on smoothed results')
    end
    % GIRF conditional on smoothed states in t_0 and shocks in t_0+1
    shocks_base= [oo_.occbin.smoother.etahat(:,t_0+1)'; zeros(options_.occbin.simul.periods,M_.exo_nbr)];
    options_.occbin.simul.SHOCKS=shocks_base;
    options_.occbin.simul.endo_init = oo_.occbin.smoother.alphahat(oo_.dr.inv_order_var,t_0);
end
options_.occbin.simul.SHOCKS=shocks_base;

[~, out_base] = occbin.solver(M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
if out_base.error_flag
    error('occbin.irfs: could not compute the solution')
end

irfs.linear = struct();
irfs.piecewise = struct();

% Get indices of shocks of interest
exo_index =zeros(size(shocknames,1),1);
for i=1:length(shocknames)
    exo_index(i) = strmatch(shocknames{i},M_.exo_names,'exact');
end

% cs=get_lower_cholesky_covariance(M_.Sigma_e,options_.add_tiny_number_to_cholesky);
% irf_shocks_indx = getIrfShocksIndx(M_, options_);

% Set shock size
if isempty(shocksize)
    shocksize = sqrt(diag(M_.Sigma_e(exo_index,exo_index)));
    if any(shocksize < 1.e-9)
        shocksize(shocksize < 1.e-9) = 0.01;
    end
end

if numel(shocksize)==1
    shocksize=repmat(shocksize,[length(shocknames),1]);
end

% Run IRFs
for sign_iter=1:length(shocksigns)
    for IRF_counter = 1:length(exo_index)
        jexo = exo_index(IRF_counter);
        if ~options_.noprint && options_.debug
            fprintf('occbin.irf: Producing GIRFs for shock %s. Simulation %d out of %d. \n',M_.exo_names{jexo},IRF_counter,size(exo_index,1));
        end
        shocks1=shocks_base;
        if ismember('pos',shocksigns{sign_iter})
            shocks1(1,jexo)=shocks_base(1,jexo)+shocksize(IRF_counter);
        elseif ismember('neg',shocksigns{sign_iter})
            shocks1(1,jexo)=shocks_base(1,jexo)-shocksize(IRF_counter);
        end
        options_.occbin.simul.SHOCKS=shocks1;
        if t_0 == 0
            options_.occbin.simul.endo_init = [];
        else
            options_.occbin.simul.endo_init = oo_.occbin.smoother.alphahat(oo_.dr.inv_order_var,t_0);
        end
        [~, out_sim] = occbin.solver(M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
        if out_sim.error_flag
            warning('occbin.irfs: simulation failed')
            skip
        end
        % Substract inital conditions + other shocks
        zdiff.linear.(shocksigns{sign_iter}) = out_sim.linear-out_base.linear;
        zdiff.piecewise.(shocksigns{sign_iter}) = out_sim.piecewise-out_base.piecewise;

        for j_endo=1:M_.endo_nbr
            if ismember('pos',shocksigns)
                irfs.piecewise.([M_.endo_names{j_endo} '_' M_.exo_names{jexo} '_' shocksigns{sign_iter}])   = zdiff.piecewise.(shocksigns{sign_iter})(:,j_endo);
                irfs.linear.([M_.endo_names{j_endo} '_' M_.exo_names{jexo} '_' shocksigns{sign_iter}])   = zdiff.linear.(shocksigns{sign_iter})(:,j_endo);
            end
        end
    end
end