function [ data, SS_out, error_flag] = solve_two_constraints(M_,dr, opts_simul_, solve_DM)
% function [ data, SS_out, error_flag] = solve_two_constraints(M_,dr, opts_simul_, solve_DM)
%
% INPUT: 
% - M_                  [structure]     Matlab's structure describing the model (M_).
% - dr                  [structure]     decision rules for the model
% - opts_simul          [structure]     Matlab's structure containing the Occbin options (opts_simul).
% - solve_DM            [double]        indicator on whether to recompute decision rules
%
% OUTPUT:
% - data                [structure]     simulation result containing fields:
%                                           - linear: paths for endogenous variables ignoring OBC (linear solution)
%                                           - piecewise: paths for endogenous variables satisfying the OBC (occbin/piecewise solution)
%                                           - ys: vector of steady state values
%                                           - regime_history: information on number and time of regime transitions
% - SS_out              [structure]     State space solution
%                                           - T: [n_vars by n_vars by n_shock_period] array of transition matrices
%                                           - R: [n_vars by n_exo by n_shock_period] array of shock response matrices
%                                           - C: [n_vars by n_shock_period] array of constants
% - error_flag          [integer]       1 if a problem was encoutered, 0 otherwise

% Original authors: Luca Guerrieri and Matteo Iacoviello 
% Original file downloaded from:
% https://www.matteoiacoviello.com/research_files/occbin_20140630.zip
% Adapted for Dynare by Dynare Team.
%
% This code is in the public domain and may be used freely.
% However the authors would appreciate acknowledgement of the source by
% citation of any of the following papers:
%
% Luca Guerrieri and Matteo Iacoviello (2015): "OccBin: A toolkit for solving
% dynamic models with occasionally binding constraints easily"
% Journal of Monetary Economics 70, 22-38

persistent DM

if isempty(DM)
    solve_DM=true;
end

data.shocks_sequence= opts_simul_.SHOCKS;                    % sequence of unforeseen shocks under which one wants to solve the model
n_periods = opts_simul_.periods;                        % simulation horizon (can be longer than the sequence of shocks defined in shockssequence; must be long enough to ensure convergence back to the reference model at the end of the simulation horizon and may need to be varied depending on the sequence of shocks).
curb_retrench = opts_simul_.curb_retrench;              % 0: updates guess based on previous iteration; 1: updates similar to Gauss-Jacobi scheme, slowing iterations down by updating guess only one period at a time
max_iter = opts_simul_.maxit;                           % maximum number of iterations allowed for the solution algorithm
endo_init = opts_simul_.endo_init;                      % initial condition for state variables, in deviation from steady state in declaration order
binding_indicator = opts_simul_.init_binding_indicator;    % initial guess for constraint violations
regime_history_guess = opts_simul_.init_regime;              % initial guess for constraint violations
periodic_solution = opts_simul_.periodic_solution;
data.exo_pos = opts_simul_.exo_pos;

n_shocks_periods = size(data.shocks_sequence,1);

if n_periods < n_shocks_periods
    n_periods = n_shocks_periods;
end
nperiods_0 = max(opts_simul_.check_ahead_periods,n_periods-n_shocks_periods);

error_flag=0;

M00_ = M_;

% ensure that all models have the same parameters
% use the parameters for the base model.

%keep the correct auxiliary regime specific parameter values

data.ys = dr.ys;

if solve_DM %recompute solution matrices
    [DM.Cbarmat ,DM.Bbarmat, DM.Abarmat, DM.Jbarmat] = occbin.get_deriv(M00_,data.ys);
    
    M10_ = M00_;
    M10_.params(M_.occbin.pswitch(1))= 1;
    [DM.Cbarmat10, DM.Bbarmat10, DM.Abarmat10, DM.Jbarmat10, DM.Dbarmat10] = occbin.get_deriv(M10_,data.ys);
    
    M01_ = M00_;
    M01_.params(M_.occbin.pswitch(2))= 1;
    [DM.Cbarmat01, DM.Bbarmat01, DM.Abarmat01, DM.Jbarmat01, DM.Dbarmat01] = occbin.get_deriv(M01_,data.ys);
    
    M11_ = M00_;
    M11_.params(M_.occbin.pswitch(1))= 1;
    M11_.params(M_.occbin.pswitch(2))= 1;
    [DM.Cbarmat11, DM.Bbarmat11, DM.Abarmat11, DM.Jbarmat11, DM.Dbarmat11] = occbin.get_deriv(M11_,data.ys);
    
    [DM.decrulea,DM.decruleb]=occbin.get_pq(dr);      
    update_flag=true;
    DM.n_vars = M00_.endo_nbr;
    DM.n_exo = M00_.exo_nbr;
else
    update_flag=false;
end

endo_names = M00_.endo_names;
exo_names =  M00_.exo_names;

init_orig_ = endo_init;

zdatapiecewise_ = zeros(n_periods,DM.n_vars);

if ~exist('binding_indicator','var')
    binding_indicator = false(nperiods_0+1,2);  % This sets the first guess for when
    % the constraints are going to hold.
    % The variable is a boolean with two columns. The first column refers to
    % constrain1_; the second to constrain2_.
    % Each row is a period in time.
    % If the boolean is true it indicates the relevant constraint is expected
    % to evaluate to true. 
    % The default initial guess is consistent with the base model always
    % holding -- equivalent to the linear solution.
else
    if size(binding_indicator,1)<(nperiods_0+1)
        binding_indicator = [binding_indicator; false(nperiods_0+1-size(binding_indicator,1),2)];
    end
end
SS_out.T = NaN(DM.n_vars,DM.n_vars,n_shocks_periods);
SS_out.R = NaN(DM.n_vars,DM.n_exo,n_shocks_periods);
SS_out.C = NaN(DM.n_vars,n_shocks_periods);
if ~exist('regime_history_','var') || isempty(regime_history_guess)
    regime_history = struct();
    guess_history = false;
else
    guess_history = true; %previous information exists
    regime_history = regime_history_guess;
end

if opts_simul_.waitbar
    hh = dyn_waitbar(0,'Occbin: Solving the model');
    set(hh,'Name','Occbin: Solving the model.');
end

for shock_period = 1:n_shocks_periods
    if opts_simul_.waitbar
        dyn_waitbar(shock_period/n_shocks_periods, hh, sprintf('Period %u of %u', shock_period,n_shocks_periods));
    end
    regime_change_this_iteration=true;
    iter = 0;
    guess_history_it = false;
    if guess_history && (shock_period<=length(regime_history_guess)) %beyond guess regime history
        guess_history_it = true;
    end
    is_periodic=false;
    binding_indicator_history={};
    max_err = NaN(max_iter,1);
    
    while (regime_change_this_iteration && iter<max_iter && ~is_periodic)
        iter = iter +1;
        if any(binding_indicator(end,:)) && nperiods_0<opts_simul_.max_periods
            binding_indicator = [binding_indicator; false(1,2)];
            nperiods_0 = nperiods_0 + 1;
            disp_verbose(['nperiods has been endogenously increased up to ' int2str(nperiods_0) '.'],opts_simul_.debug)
        end
        if size(binding_indicator,1)<(nperiods_0 + 1)
            binding_indicator=[binding_indicator; false(nperiods_0 + 1-size(binding_indicator,1),2)];
        end
        binding_indicator_history{iter}=binding_indicator;
        
        if iter==1 && guess_history_it
            regime_1 = regime_history_guess(shock_period).regime1;
            regime_start_1 = regime_history_guess(shock_period).regimestart1;
            binding_indicator(:,1) = regime_1(end);
            for ir=1:length(regime_1)-1
                binding_indicator(regime_start_1(ir):regime_start_1(ir+1)-1,1) = regime_1(ir);
            end
            regime_2 = regime_history_guess(shock_period).regime2;
            regime_start_2 = regime_history_guess(shock_period).regimestart2;
            binding_indicator(:,2) = regime_2(end);
            for ir=1:length(regime_2)-1
                binding_indicator(regime_start_2(ir):regime_start_2(ir+1)-1,2) = regime_2(ir);
            end
            nperiods_0 = size(binding_indicator,1)-1; %if history is present, update may be required
        end
        % analyse violvec and isolate contiguous periods in the other regime.
        [regime_1, regime_start_1, error_code_period(1)]=occbin.map_regime(binding_indicator(:,1),opts_simul_.debug);
        regime_history(shock_period).regime1 = regime_1;
        regime_history(shock_period).regimestart1 = regime_start_1;
        [regime_2, regime_start_2, error_code_period(2)]=occbin.map_regime(binding_indicator(:,2),opts_simul_.debug);
        regime_history(shock_period).regime2 = regime_2;
        regime_history(shock_period).regimestart2 = regime_start_2;
        if shock_period==1 || shock_period>1 && any(data.shocks_sequence(shock_period,:)) % first period or shock happening
            if iter==1 && opts_simul_.reset_regime_in_new_period
                binding_indicator=false(size(binding_indicator));
                binding_indicator_history{iter}=binding_indicator;
                % analyse violvec and isolate contiguous periods in the other regime.
                [regime_1, regime_start_1, error_code_period(1)]=occbin.map_regime(binding_indicator(:,1),opts_simul_.debug);
                regime_history(shock_period).regime1 = regime_1;
                regime_history(shock_period).regimestart1 = regime_start_1;
                [regime_2, regime_start_2, error_code_period(2)]=occbin.map_regime(binding_indicator(:,2),opts_simul_.debug);
                regime_history(shock_period).regime2 = regime_2;
                regime_history(shock_period).regimestart2 = regime_start_2;
            end
            Tmax=max([regime_start_1(end) regime_start_2(end)])-1;
            [zdatalinear_, SS_out.T(:,:,shock_period), SS_out.R(:,:,shock_period), SS_out.C(:,shock_period), SS, update_flag]=occbin.mkdatap_anticipated_2constraints_dyn(nperiods_0,...
                DM,Tmax,...
                binding_indicator,...
                data.exo_pos,data.shocks_sequence(shock_period,:),endo_init, update_flag);
            
            [binding, relax, err]=feval([M_.fname,'.occbin_difference'],zdatalinear_+repmat(dr.ys',size(zdatalinear_,1),1),M_.params,dr.ys);
            binding_constraint_new=[binding.constraint_1;binding.constraint_2];
            relaxed_constraint_new = [relax.constraint_1;relax.constraint_2];
            
            err_binding_constraint_new = [err.binding_constraint_1; err.binding_constraint_2];
            err_relaxed_constraint_new = [err.relax_constraint_1; err.relax_constraint_2];
            
            % check if changes_
            if any(binding_constraint_new & ~binding_indicator(:)) || any(relaxed_constraint_new & binding_indicator(:))
                err_violation = err_binding_constraint_new(binding_constraint_new & ~binding_indicator(:));
                err_relax = err_relaxed_constraint_new(relaxed_constraint_new & binding_indicator(:));
                max_err(iter) = max(abs([err_violation;err_relax]));
                regime_change_this_iteration = true;
            else
                regime_change_this_iteration = false;
                max_err(iter) = 0;
            end
            
            if curb_retrench   % apply Gauss-Seidel idea of slowing down the change in the guess
                % for the constraint -- only relax one
                % period at a time starting from the last
                % one when each of the constraints is true.
                retrench = false(numel(binding_indicator),1);
                max_relax_constraint_1=find(relax.constraint_1 & binding_indicator(:,1),1,'last');
                if ~isempty(max_relax_constraint_1) && find(relax.constraint_1,1,'last')>=find(binding_indicator(:,1),1,'last')
                    retrench(max_relax_constraint_1) = true;
                end
                max_relax_constraint_2=find(relax.constraint_2 & binding_indicator(:,2),1,'last');
                if ~isempty(max_relax_constraint_2) && find(relax.constraint_2,1,'last')>=find(binding_indicator(:,2),1,'last')
                    retrench(max_relax_constraint_2+nperiods_0+1) = true;
                end
                binding_indicator = (binding_indicator(:) | binding_constraint_new) & ~ retrench;
            else
                binding_indicator= (binding_indicator(:) | binding_constraint_new) & ~(binding_indicator(:) & relaxed_constraint_new);
            end
            binding_indicator = reshape(binding_indicator,nperiods_0+1,2);
            
            if iter>1 && regime_change_this_iteration
                is_periodic=false(1,iter-1);
                for kiter=1:iter-1
                    vvv = [binding_indicator_history{kiter}; false(size(binding_indicator,1)- size(binding_indicator_history{kiter},1), 2)];
                    is_periodic(kiter) = isequal(vvv, binding_indicator);
                end
                is_periodic_all = is_periodic;
                is_periodic =  any(is_periodic);
                if is_periodic && periodic_solution
                    [min_err,index_min_err]=min(max_err(find(is_periodic_all,1):end));
                    inx = find(is_periodic_all,1):iter;
                    inx = inx(index_min_err);
                    binding_indicator=binding_indicator_history{inx}; %select regime history with same result, but smallest error
                    if inx<iter %update if needed
                        [regime_1, regime_start_1, error_code_period(1)]=occbin.map_regime(binding_indicator(:,1),opts_simul_.debug);
                        regime_history(shock_period).regime1 = regime_1;
                        regime_history(shock_period).regimestart1 = regime_start_1;
                        [regime_2, regime_start_2, error_code_period(2)]=occbin.map_regime(binding_indicator(:,2),opts_simul_.debug);
                        regime_history(shock_period).regime2 = regime_2;
                        regime_history(shock_period).regimestart2 = regime_start_2;
                        Tmax=max([regime_start_1(end) regime_start_2(end)])-1;
                        % get the hypothesized piece wise linear solution
                        [zdatalinear_, SS_out.T(:,:,shock_period), SS_out.R(:,:,shock_period), SS_out.C(:,shock_period), SS, update_flag]=occbin.mkdatap_anticipated_2constraints_dyn(nperiods_0,DM,...
                            Tmax,...
                            binding_indicator,...
                            data.exo_pos,data.shocks_sequence(shock_period,:),endo_init,update_flag);                        
                    end
                end
                
            end
        else
            regime_change_this_iteration= false;
            zdatalinear_(1:end-1,:)=zdatalinear_(2:end,:);
            zdatalinear_(end,:) = DM.decrulea*zdatalinear_(end-1,:)';
            if length(SS)>1
                SS=SS(2:end);
            else
                SS=[];
            end
            if isempty(SS)
                SS_out.T(:,:,shock_period)= DM.decrulea;
                SS_out.R(:,:,shock_period)= DM.decruleb;
                SS_out.C(:,shock_period)= 0;
            else
                SS_out.T(:,:,shock_period)= SS(1).T;
                SS_out.R(:,:,shock_period)= SS(1).R;
                SS_out.C(:,shock_period)= SS(1).C;
            end
            binding_indicator_history{iter}=binding_indicator;
        end
    end
    if regime_change_this_iteration
        if max_iter>opts_simul_.algo_truncation
            disp_verbose(['occbin solver: period ' int2str(shock_period) ':'],opts_simul_.debug)
            if is_periodic
                disp_verbose('Occbin solver loops between two regimes.',opts_simul_.debug)
                if periodic_solution
                    disp_verbose(['Max error:' num2str(min_err) '.'],opts_simul_.debug)
                else
                    error_flag = 310;
                    if opts_simul_.waitbar; dyn_waitbar_close(hh); end
                    return;
                end
            else
                disp_verbose('Did not converge -- increase maxit.',opts_simul_.debug)
                error_flag = 311;
                if opts_simul_.waitbar; dyn_waitbar_close(hh); end
                return;
            end
        else
            % if max_iter <= truncation, we force indicator to equal the
            % last guess
            binding_indicator = binding_indicator_history{end};
        end
    end
    if any(error_code_period)
        disp_verbose('Increase nperiods.',opts_simul_.debug)
        error_flag = 312;
        if opts_simul_.waitbar; dyn_waitbar_close(hh); end
        return;
    end
    
    endo_init = zdatalinear_(1,:);
    zdatapiecewise_(shock_period,:)=endo_init;
    endo_init= endo_init';
    
    % update the guess for constraint violations for next period
    % update is consistent with expecting no additional shocks next period
    binding_indicator=[binding_indicator(2:end,:); false(1,2)];
    
end

zdatapiecewise_(shock_period+1:end,:)=zdatalinear_(2:n_periods-shock_period+1,:);

data.piecewise=zdatapiecewise_;
data.regime_history=regime_history;

if ~opts_simul_.piecewise_only
    % get the linear responses
    data.linear = occbin.mkdata(n_periods,DM.decrulea,DM.decruleb,endo_names,exo_names,[],data.exo_pos,data.shocks_sequence,init_orig_);
end

if opts_simul_.waitbar
    dyn_waitbar_close(hh); 
end
