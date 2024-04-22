function oo_=disp_moments(y,var_list,M_,options_,oo_)
% function disp_moments(y,var_list,M_,options_,oo_)
% Displays moments of simulated variables
% INPUTS
%   y                   [double]       nvar*nperiods vector of simulated variables.
%   var_list            [char]         nvar character array with names of variables.
%   M_                  [structure]    Dynare's model structure
%   oo_                 [structure]    Dynare's results structure
%   options_            [structure]    Dynare's options structure
%
% OUTPUTS
%   oo_                 [structure]    Dynare's results structure,

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

zero_moments_tolerance=1e-10;
warning_old_state = warning;
warning off

if isempty(var_list)
    var_list = M_.endo_names(1:M_.orig_endo_nbr);
end

nvar = length(var_list);
ivar=zeros(nvar,1);
for i=1:nvar
    i_tmp = strmatch(var_list{i}, M_.endo_names, 'exact');
    if isempty(i_tmp)
       error('The variable %s specified is not an endogenous variable',var_list{i});
    else
        ivar(i) = i_tmp;
    end
end

y = y(ivar,options_.drop+1:end)';

ME_present=0;
if ~all(M_.H==0)
    if isoctave && octave_ver_less_than('8.4') %Octave bug #60347
        [observable_pos_requested_vars, index_subset, index_observables] = intersect_stable(ivar, options_.varobs_id);
    else
        [observable_pos_requested_vars, index_subset, index_observables] = intersect(ivar, options_.varobs_id, 'stable');
    end
    if ~isempty(observable_pos_requested_vars)
        ME_present=1;
        i_ME = setdiff(1:size(M_.H,1),find(diag(M_.H) == 0)); % find ME with 0 variance
        chol_S = chol(M_.H(i_ME,i_ME)); %decompose rest
        shock_mat=zeros(options_.periods,size(M_.H,1)); %initialize
        shock_mat(:,i_ME)=randn(length(i_ME),options_.periods)'*chol_S;
        y_ME = y(:,index_subset)+shock_mat(options_.drop+1:end,index_observables);
        y_ME_only = shock_mat(options_.drop+1:end,index_observables);
        m_ME = mean(y_ME);
        y_ME=get_filtered_time_series(y_ME,m_ME,options_);
        y_ME_only_filtered=get_filtered_time_series(y_ME_only,mean(y_ME_only),options_);
        s2_ME = mean(y_ME.*y_ME);
        s_ME = sqrt(s2_ME);
        zero_variance_ME_var_index=index_subset(abs(s_ME')<zero_moments_tolerance);
    end
end


m = mean(y);

% filter series
y=get_filtered_time_series(y,m,options_);

s2 = mean(y.*y);
s = sqrt(s2);
oo_.mean = transpose(m);
oo_.var = y'*y/size(y,1);
oo_.skewness = (mean(y.^3)./s2.^1.5)';
oo_.kurtosis = (mean(y.^4)./(s2.*s2)-3)';

zero_variance_var_index=find(abs(s)<zero_moments_tolerance);
oo_.skewness(zero_variance_var_index)=NaN;
oo_.kurtosis(zero_variance_var_index)=NaN;
s(zero_variance_var_index)=0;
s2(zero_variance_var_index)=0;
oo_.var(zero_variance_var_index,:)=0;
oo_.var(:,zero_variance_var_index)=0;


labels=get_labels_transformed_vars(M_.endo_names,ivar,options_,false);
labels_TeX=get_labels_transformed_vars(M_.endo_names_tex,ivar,options_,true);

if ~options_.nomoments
    z = [ m' s' s2' oo_.skewness oo_.kurtosis ];
    title='MOMENTS OF SIMULATED VARIABLES';
    title=add_filter_subtitle(title, options_);
    headers = {'VARIABLE'; 'MEAN'; 'STD. DEV.'; 'VARIANCE'; 'SKEWNESS'; 'KURTOSIS'};
    dyntable(options_, title, headers, labels, z, cellofchararraymaxlength(labels)+2, 16, 6);
    if options_.TeX
        dyn_latex_table(M_, options_, title, 'sim_moments', headers, labels_TeX, z, cellofchararraymaxlength(labels)+2, 16, 6);
    end
end

if ~options_.nocorr
    corr = (y'*y/size(y,1))./(s'*s);
    corr(zero_variance_var_index,:)=NaN;
    corr(:,zero_variance_var_index)=NaN;
    if options_.contemporaneous_correlation
        oo_.contemporaneous_correlation = corr;
    end
    if ~options_.noprint
        title = 'CORRELATION OF SIMULATED VARIABLES';
        title=add_filter_subtitle(title,options_);
        headers = vertcat('VARIABLE', M_.endo_names(ivar));
        dyntable(options_, title, headers, labels, corr, cellofchararraymaxlength(labels)+2, 8, 4);
        if options_.TeX
            headers = vertcat('VARIABLE', M_.endo_names_tex(ivar));
            lh = cellofchararraymaxlength(labels)+2;
            dyn_latex_table(M_, options_, title, 'sim_corr_matrix', headers, labels_TeX, corr, lh, 8,4);
        end
    end
end

if ~options_.noprint && length(options_.conditional_variance_decomposition)
    fprintf('\nSTOCH_SIMUL: conditional_variance_decomposition requires theoretical moments, i.e. periods=0.\n')
end

ar = options_.ar;
if ar > 0
    autocorr = [];
    for i=1:ar
        oo_.autocorr{i} = y(ar+1:end,:)'*y(ar+1-i:end-i,:)./((size(y,1)-ar)*std(y(ar+1:end,:))'*std(y(ar+1-i:end-i,:)));
        oo_.autocorr{i}(zero_variance_var_index,:)=NaN;
        oo_.autocorr{i}(:,zero_variance_var_index)=NaN;
        autocorr = [ autocorr diag(oo_.autocorr{i}) ];
    end
    if ~options_.noprint
        title = 'AUTOCORRELATION OF SIMULATED VARIABLES';
        title=add_filter_subtitle(title,options_);
        headers = vertcat('VARIABLE', cellstr(int2str([1:ar]')));
        dyntable(options_, title, headers, labels, autocorr, cellofchararraymaxlength(labels)+2, 8, 4);
        if options_.TeX
            headers = vertcat('VARIABLE', cellstr(int2str([1:ar]')));
            lh = cellofchararraymaxlength(labels)+2;
            dyn_latex_table(M_, options_, title, 'sim_autocorr_matrix', headers, labels_TeX, autocorr, cellofchararraymaxlength(labels_TeX)+2, 8, 4);
        end
    end

end


if ~options_.nodecomposition
    if M_.exo_nbr == 1
        oo_.variance_decomposition = 100*ones(nvar,1);
    else
        oo_.variance_decomposition=zeros(nvar,M_.exo_nbr);
        %get starting values
        if isempty(M_.endo_histval)
            y0 = oo_.dr.ys;
        else
            if options_.loglinear
                y0 = log_variable(1:M_.endo_nbr,M_.endo_histval,M_);
            else
                y0 = M_.endo_histval;
            end
        end
        %back out shock matrix used for generating y
        i_exo_var = setdiff(1:M_.exo_nbr,find(diag(M_.Sigma_e) == 0)); % find shocks with 0 variance
        chol_S = chol(M_.Sigma_e(i_exo_var,i_exo_var)); %decompose rest
        shock_mat=zeros(options_.periods,M_.exo_nbr); %initialize
        shock_mat(:,i_exo_var)=oo_.exo_simul(:,i_exo_var)/chol_S; %invert construction of oo_.exo_simul from simult.m

        for shock_iter=1:length(i_exo_var)
            temp_shock_mat=zeros(size(shock_mat));
            temp_shock_mat(:,i_exo_var(shock_iter))=shock_mat(:,i_exo_var(shock_iter));
            temp_shock_mat(:,i_exo_var) = temp_shock_mat(:,i_exo_var)*chol_S;
            y_sim_one_shock = simult_(M_,options_,y0,oo_.dr,temp_shock_mat,options_.order);
            y_sim_one_shock=y_sim_one_shock(ivar,1+options_.drop+1:end)';
            y_sim_one_shock=get_filtered_time_series(y_sim_one_shock,mean(y_sim_one_shock),options_);
            oo_.variance_decomposition(:,i_exo_var(shock_iter))=var(y_sim_one_shock)./s2*100;            
        end
        oo_.variance_decomposition(zero_variance_var_index,:)=NaN;
        if ME_present
            oo_.variance_decomposition_ME=oo_.variance_decomposition(index_subset,:)...
                .*repmat((s2(index_subset)./s2_ME)',1,M_.exo_nbr);
            oo_.variance_decomposition_ME(:,end+1)=var(y_ME_only_filtered)./s2_ME*100;
            oo_.variance_decomposition_ME(ismember(observable_pos_requested_vars,intersect(zero_variance_ME_var_index,zero_variance_var_index)),:)=NaN;
            oo_.variance_decomposition_ME(ismember(observable_pos_requested_vars,setdiff(zero_variance_var_index,zero_variance_ME_var_index)),1:end-1)=0;
            oo_.variance_decomposition_ME(ismember(observable_pos_requested_vars,setdiff(zero_variance_var_index,zero_variance_ME_var_index)),end)=1;
        end
        if ~options_.noprint %options_.nomoments == 0
            skipline()
            title='VARIANCE DECOMPOSITION SIMULATING ONE SHOCK AT A TIME (in percent)';
            title=add_filter_subtitle(title,options_);
            headers = M_.exo_names;
            headers = vertcat(' ', headers);
            lh = cellofchararraymaxlength(M_.endo_names(ivar))+2;
            dyntable(options_, title, vertcat(headers, 'Tot. lin. contr.'), ...
                     labels, [oo_.variance_decomposition sum(oo_.variance_decomposition,2)], lh, 8, 2);
            if ME_present
                headers_ME = vertcat(headers, 'ME');
                labels_obs=get_labels_transformed_vars(M_.endo_names,ivar(index_subset),options_,false);
                dyntable(options_, [title,' WITH MEASUREMENT ERROR'], vertcat(headers_ME, 'Tot. lin. contr.'), labels_obs, ...
                         [oo_.variance_decomposition_ME sum(oo_.variance_decomposition_ME, 2)], lh, 8, 2);
            end
            if options_.TeX
                headers = M_.exo_names_tex;
                headers = vertcat(' ', headers);
                labels = M_.endo_names_tex(ivar);
                lh = cellofchararraymaxlength(labels)+2;
                dyn_latex_table(M_, options_, title, 'sim_var_decomp', vertcat(headers, 'Tot. lin. contr.'), ...
                                labels_TeX, [oo_.variance_decomposition sum(oo_.variance_decomposition, 2)], lh, 8, 2);
                if ME_present
                    headers_ME = vertcat(headers, 'ME');
                    labels_obs_TeX=get_labels_transformed_vars(M_.endo_names_tex,ivar(index_subset),options_,true);
                    dyn_latex_table(M_, options_, [title, ' WITH MEASUREMENT ERROR'], 'sim_var_decomp_ME', ...
                                    vertcat(headers_ME, 'Tot. lin. contr.'), ...
                                    labels_obs_TeX, ...
                                    [oo_.variance_decomposition_ME sum(oo_.variance_decomposition_ME, 2)], lh, 8, 2);
                end
            end

            if options_.order == 1
                fprintf('Note: numbers do not add up to 100 due to non-zero correlation of simulated shocks in small samples\n\n')
            else
                fprintf('Note: numbers do not add up to 100 due to i) non-zero correlation of simulated shocks in small samples and ii) nonlinearity\n\n')
            end
        end

    end
end

warning(warning_old_state);
end

function y = get_filtered_time_series(y, m, options_)

if options_.hp_filter && ~options_.one_sided_hp_filter  && ~options_.bandpass.indicator
    [~,y] = sample_hp_filter(y,options_.hp_filter);
elseif ~options_.hp_filter && options_.one_sided_hp_filter && ~options_.bandpass.indicator
    [~,y] = one_sided_hp_filter(y,options_.one_sided_hp_filter);
elseif ~options_.hp_filter && ~options_.one_sided_hp_filter && options_.bandpass.indicator
    data_temp=dseries(y,'0q1');
    data_temp=baxter_king_filter(data_temp,options_.bandpass.passband(1),options_.bandpass.passband(2),options_.bandpass.K);
    y=data_temp.data;
elseif ~options_.hp_filter && ~options_.one_sided_hp_filter  && ~options_.bandpass.indicator
    y = bsxfun(@minus, y, m);
else
    error('disp_moments:: You cannot use more than one filter at the same time')
end

end

