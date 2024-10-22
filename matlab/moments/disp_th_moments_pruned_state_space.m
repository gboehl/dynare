function oo_=disp_th_moments_pruned_state_space(dr,M_,options_,i_var,oo_)
% oo_=disp_th_moments_pruned_state_space(dr,M_,options_,i_var,oo_)
% Display theoretical moments of variables based on (second or third order)
% pruned state-space
%
% INPUTS:
% dr :          [struct]    Dynare decision rules structure
% M_            [struct]    structure describing the Model
% options_      [struct]    structure describing the options
% i_var         [double]    Index of requested variables in declaration order
% oo_           [struct]    structure describing the Model
%
% OUTPUTS: 
%           gamma_y                                 [cell]      Matlab cell of nar+1 arrays, where nar is the order of the autocorrelation function.
%           gamma_y{1}                              [double]    Covariance matrix.
%           gamma_y{i+1}                            [double]    Autocorrelation function (for i=1,...,options_.ar).
%           mean                                    [vector]    Unconditional mean
%           var                                     [matrix]    Unconditional covariance matrix
%           autocorr                                [cell]      Cell storing the theoretical autocorrelation
%           contemporaneous_correlation             [matrix]    matrix of contemporaneous correlations
%           autocorr                                [cell]      Cell storing the theoretical autocorrelation
%           variance_decomposition                  [matrix]    Unconditional variance decomposition matrix
%           variance_decomposition_ME               [matrix]    Unconditional variance decomposition matrix with measurement error
%           conditional_variance_decomposition      [array]     Conditional variance decomposition array
%           conditional_variance_decomposition_ME   [array]     Conditional variance decomposition array with measurement error

% Copyright © 2020-2023 Dynare Team
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


if options_.one_sided_hp_filter || options_.hp_filter || options_.bandpass.indicator
    error('disp_th_moments:: theoretical moments incompatible with filtering. Use simulated moments instead')
end

nvars=length(i_var);
obs_var=NaN(nvars,1);
for i=1:nvars
    obs_var(i,1) = find(strcmp(M_.endo_names(i_var(i),:), M_.endo_names(dr.order_var)));
end

pruned_state_space = pruned_SS.pruned_state_space_system(M_, options_, dr, obs_var, options_.ar, 1, 0);


m = pruned_state_space.E_y;

oo_.gamma_y{1} = pruned_state_space.Var_y;

i1 = find(abs(diag(oo_.gamma_y{1})) > 1e-12);
s2 = diag(oo_.gamma_y{1});
sd = sqrt(s2);

z = [ m sd s2 ];
oo_.mean = m;
oo_.var = oo_.gamma_y{1};

if ~options_.noprint %options_.nomoments == 0
    title='THEORETICAL MOMENTS BASED ON PRUNED STATE SPACE';
    headers={'VARIABLE','MEAN','STD. DEV.','VARIANCE'};
    labels = M_.endo_names(i_var,:);
    lh = cellofchararraymaxlength(labels)+2;
    dyntable(options_,title,headers,labels,z,lh,11,4);
    if options_.TeX
        labels = M_.endo_names_tex(i_var,:);
        lh = cellofchararraymaxlength(labels)+2;
        dyn_latex_table(M_,options_,title,'th_moments',headers,labels,z,lh,11,4);
    end        
end

if isempty(i1)
    disp_verbose(' ',~options_.noprint)
    disp_verbose('All endogenous are constant or non stationary, not displaying correlations and auto-correlations',~options_.noprint)
    disp_verbose(' ',~options_.noprint)
    return;
end

if options_.nocorr == 0 % && size(stationary_vars, 1) > 0
    corr=pruned_state_space.Corr_y;
    if options_.contemporaneous_correlation 
        oo_.contemporaneous_correlation = corr;
    end
    if ~options_.noprint
        skipline()
        title='MATRIX OF CORRELATIONS BASED ON PRUNED STATE SPACE';            
        labels = M_.endo_names(i_var,:);
        headers = ['Variables';labels];
        lh = cellofchararraymaxlength(labels)+2;
        dyntable(options_,title,headers,labels,corr,lh,8,4);
        if options_.TeX
            labels = M_.endo_names_tex(i_var,:);
            headers=['Variables';labels];
            lh = cellofchararraymaxlength(labels)+2;
            dyn_latex_table(M_,options_,title,'th_corr_matrix',headers,labels,corr,lh,8,4);
        end
    end
end
if options_.ar > 0 %&& size(stationary_vars, 1) > 0
    z=NaN(length(i1),options_.ar);
    for i=1:options_.ar
        oo_.gamma_y{i+1} = pruned_state_space.Corr_yi(:,:,i);
        oo_.autocorr{i} = oo_.gamma_y{i+1};
        z(:,i) = diag(oo_.gamma_y{i+1}(i1,i1));
    end
    if ~options_.noprint    
        skipline()    
        title='COEFFICIENTS OF AUTOCORRELATION BASED ON PRUNED STATE SPACE';            
        labels = M_.endo_names(i_var(i1),:);
        headers = ['Order ';cellstr(int2str([1:options_.ar]'))];
        lh = cellofchararraymaxlength(labels)+2;
        dyntable(options_,title,headers,labels,z,lh,8,4);
        if options_.TeX
            labels = M_.endo_names_tex(i_var(i1),:);
            lh = cellofchararraymaxlength(labels)+2;
            dyn_latex_table(M_,options_,title,'th_autocorr_matrix',headers,labels,z,lh,8,4);
        end
    end  
end

if options_.order==2 && ~options_.nodecomposition && M_.exo_nbr > 1% do variance decomposition
    index_stationary_vars = dr.inv_order_var(i_var); %no nonstationary variables in pruning
    ghu_states_only = zeros(M_.nspred,M_.exo_nbr);
    ghu_states_only(1:M_.nspred,:) = dr.ghu(M_.nstatic+(1:M_.nspred),:); %get shock impact on states only
    [A] = kalman_transition_matrix(dr,M_.nstatic+(1:M_.nspred)',1:M_.nspred);
    oo_.gamma_y{options_.ar+2}=compute_variance_decomposition(M_,options_,s2,A,dr.ghx(index_stationary_vars,:),dr.ghu,ghu_states_only,1:length(i_var),index_stationary_vars,nvars);

    [ME_present,observable_pos_requested_vars,index_subset,index_observables]=check_measurement_error_requested_vars(M_,options_,i_var);
    %store unconditional variance decomposition
    oo_.variance_decomposition=100*oo_.gamma_y{options_.ar+2};
    if ME_present
        ME_Variance=diag(M_.H);
        oo_.variance_decomposition_ME=oo_.variance_decomposition(index_subset,:).*repmat(diag(oo_.var(index_subset,index_subset))./(diag(oo_.var(index_subset,index_subset))+ME_Variance(index_observables)),1,M_.exo_nbr);
        oo_.variance_decomposition_ME(:,end+1)=100-sum(oo_.variance_decomposition_ME,2);
    end
    if ~options_.noprint %options_.nomoments == 0
        display_unconditional_variance_decomposition(M_,options_,oo_,i_var,1:length(i_var),index_subset,ME_present)
    end

    %% Conditional variance decomposition
    conditional_variance_steps = options_.conditional_variance_decomposition;
    if ~isempty(conditional_variance_steps)
        [oo_.conditional_variance_decomposition, oo_.conditional_variance_decomposition_ME] = ...
            conditional_variance_decomposition(M_,options_,dr, conditional_variance_steps, i_var);
        if ~options_.noprint
            display_conditional_variance_decomposition(oo_.conditional_variance_decomposition, conditional_variance_steps, i_var, M_, options_);
            if ME_present
                display_conditional_variance_decomposition(oo_.conditional_variance_decomposition_ME, conditional_variance_steps, ...
                    observable_pos_requested_vars, M_, options_);
            end
        end
    end
end