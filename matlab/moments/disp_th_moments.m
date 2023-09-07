function oo_ = disp_th_moments(dr, var_list, M_, options_, oo_)
%oo_ = disp_th_moments(dr, var_list, M_, options_, oo_)
% Display theoretical moments of variables
% INPUTS:
% dr :          [struct]    Dynare decision rules structure
% var_list      [cell]      list of variables considered
% M_            [struct]    structure describing the Model
% options_      [struct]    structure describing the options
% oo_           [struct]    structure describing the Model
%
% OUTPUTS: 
% oo_           [struct]    structure describing the Model, containing       
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

% Copyright © 2001-2023 Dynare Team
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

nodecomposition = options_.nodecomposition;
if options_.one_sided_hp_filter
    error('disp_th_moments:: theoretical moments incompatible with one-sided HP filter. Use simulated moments instead')
end
if isempty(var_list)
    var_list = M_.endo_names(1:M_.orig_endo_nbr);
end
nvar = length(var_list);
ivar=zeros(nvar,1);
for i=1:nvar
    i_tmp = strmatch(var_list{i}, M_.endo_names, 'exact');
    if isempty(i_tmp)
        error ('One of the variable specified does not exist');
    else
        ivar(i) = i_tmp;
    end
end

[oo_.gamma_y,stationary_vars] = th_autocovariances(dr, ivar, M_, options_, nodecomposition);
m = dr.ys(ivar);
non_stationary_vars = setdiff(1:length(ivar),stationary_vars);
m(non_stationary_vars) = NaN;

i1 = find(abs(diag(oo_.gamma_y{1})) > 1e-12);
s2 = diag(oo_.gamma_y{1});
sd = sqrt(s2);
if options_.order == 2 && ~M_.hessian_eq_zero
    m = m+oo_.gamma_y{options_.ar+3};
end

z = [ m sd s2 ];
oo_.mean = m;
oo_.var = oo_.gamma_y{1};


if size(stationary_vars, 1) > 0
    if ~options_.noprint
        if options_.order == 2
            title = 'APPROXIMATED THEORETICAL MOMENTS';
        else
            title = 'THEORETICAL MOMENTS';
        end
        title = add_filter_subtitle(title, options_);
        headers = {'VARIABLE';'MEAN';'STD. DEV.';'VARIANCE'};
        labels=get_labels_transformed_vars(M_.endo_names,ivar,options_,false);
        lh = cellofchararraymaxlength(labels)+2;
        dyntable(options_, title, headers, labels, z, lh, 11, 4);
        if options_.TeX
            labels=get_labels_transformed_vars(M_.endo_names_tex,ivar,options_,true);
            lh = cellofchararraymaxlength(labels)+2;
            dyn_latex_table(M_, options_, title, 'th_moments', headers, labels, z, lh, 11, 4);
        end
    end

    [ME_present,observable_pos_requested_vars,index_subset,index_observables]=check_measurement_error_requested_vars(M_,options_,ivar);
    %store unconditional variance decomposition
    if ~nodecomposition
        oo_.variance_decomposition=100*oo_.gamma_y{options_.ar+2};
        if ME_present
            ME_Variance=diag(M_.H);
            oo_.variance_decomposition_ME=oo_.variance_decomposition(index_subset,:).*repmat(diag(oo_.var(index_subset,index_subset))./(diag(oo_.var(index_subset,index_subset))+ME_Variance(index_observables)),1,M_.exo_nbr);
            oo_.variance_decomposition_ME(:,end+1)=100-sum(oo_.variance_decomposition_ME,2);
        end
    end
    if ~options_.noprint %options_.nomoments == 0
        if M_.exo_nbr > 1 && ~nodecomposition
            display_unconditional_variance_decomposition(M_,options_,oo_,ivar,stationary_vars,index_subset,ME_present)
        end
    end
end
%% Conditional variance decomposition
conditional_variance_steps = options_.conditional_variance_decomposition;
if ~isempty(conditional_variance_steps)
    [oo_.conditional_variance_decomposition, oo_.conditional_variance_decomposition_ME] = ...
        conditional_variance_decomposition(M_,options_,dr, conditional_variance_steps, ivar);
    if ~options_.noprint
        display_conditional_variance_decomposition(oo_.conditional_variance_decomposition, conditional_variance_steps, ivar, M_, options_);
        if ME_present
            display_conditional_variance_decomposition(oo_.conditional_variance_decomposition_ME, conditional_variance_steps, ...
                observable_pos_requested_vars, M_, options_);
        end
    end
end

if isempty(i1)
    if ~options_.noprint
        skipline()
        disp('All endogenous are constant or non stationary, not displaying correlations and auto-correlations')
        skipline()
    end
    return
end

if ~options_.nocorr && size(stationary_vars, 1)>0
    corr = NaN(size(oo_.gamma_y{1}));
    corr(i1,i1) = oo_.gamma_y{1}(i1,i1)./(sd(i1)*sd(i1)');
    if options_.contemporaneous_correlation
        oo_.contemporaneous_correlation = corr;
    end
    if ~options_.noprint
        skipline()
        if options_.order==2
            title = 'APPROXIMATED MATRIX OF CORRELATIONS';
        else
            title = 'MATRIX OF CORRELATIONS';
        end
        title = add_filter_subtitle(title, options_);
        labels=get_labels_transformed_vars(M_.endo_names,ivar(i1),options_,false);
        headers = vertcat('Variables', labels);
        lh = cellofchararraymaxlength(labels)+2;
        dyntable(options_, title, headers, labels, corr(i1,i1), lh, 8, 4);
        if options_.TeX
            labels=get_labels_transformed_vars(M_.endo_names_tex,ivar(i1),options_,true);
            headers = vertcat('Variables', labels);
            lh = cellofchararraymaxlength(labels)+2;
            dyn_latex_table(M_, options_, title, 'th_corr_matrix', headers, labels, corr(i1,i1), lh, 8, 4);
        end
    end
end

if options_.ar > 0 && size(stationary_vars, 1) > 0
    z=NaN(length(i1),options_.ar);
    for i=1:options_.ar
        oo_.autocorr{i} = oo_.gamma_y{i+1};
        z(:,i) = diag(oo_.gamma_y{i+1}(i1,i1));
    end
    if ~options_.noprint
        skipline()
        if options_.order == 2
            title = 'APPROXIMATED COEFFICIENTS OF AUTOCORRELATION';
        else
            title = 'COEFFICIENTS OF AUTOCORRELATION';
        end
        title = add_filter_subtitle(title, options_);
        labels=get_labels_transformed_vars(M_.endo_names,ivar(i1),options_,false);
        headers = vertcat('Order ', cellstr(int2str((1:options_.ar)')));
        lh = cellofchararraymaxlength(labels)+2;
        dyntable(options_, title, headers, labels, z, lh, 8, 4);
        if options_.TeX
            labels=get_labels_transformed_vars(M_.endo_names_tex,ivar(i1),options_,true);
            headers = vertcat('Order ', cellstr(int2str((1:options_.ar)')));
            lh = cellofchararraymaxlength(labels)+2;
            dyn_latex_table(M_, options_, title, 'th_autocorr_matrix', headers, labels, z, lh, 8, 4);
        end
    end
end
