function [dataMoments, m_data] = data_moments(data, oo_, matched_moments_, options_mom_)
% [dataMoments, m_data] = data_moments(data, oo_, matched_moments_, options_mom_)
% This function computes the user-selected empirical moments from data
% =========================================================================
% INPUTS
%  o data                    [T x varobs_nbr]  data set
%  o oo_:                    [structure]       storage for results
%  o matched_moments_:       [structure]       information about selected moments to match in estimation
%  o options_mom_:           [structure]       information about all settings (specified by the user, preprocessor, and taken from global options_)
% -------------------------------------------------------------------------
% OUTPUTS
%  o dataMoments             [numMom x 1]       mean of selected empirical moments
%  o m_data                  [T x numMom]       selected empirical moments at each point in time
% -------------------------------------------------------------------------
% This function is called by
%  o mom.run.m
%  o mom.objective_function.m
% =========================================================================
% Copyright (C) 2020-2021 Dynare Team
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
% -------------------------------------------------------------------------
% Author(s): 
% o Willi Mutschler (willi@mutschler.eu)
% o Johannes Pfeifer (jpfeifer@uni-koeln.de)
% =========================================================================

% Initialization
T = size(data,1); % Number of observations (T)
dataMoments = NaN(options_mom_.mom.mom_nbr,1);
m_data = NaN(T,options_mom_.mom.mom_nbr);
% Product moment for each time period, i.e. each row t contains y_t1(l1)^p1*y_t2(l2)^p2*...
% note that here we already are able to treat leads and lags and any power product moments
for jm = 1:options_mom_.mom.mom_nbr
    vars     = oo_.dr.inv_order_var(matched_moments_{jm,1})';
    leadlags = matched_moments_{jm,2}; % lags are negative numbers and leads are positive numbers
    powers   = matched_moments_{jm,3};
    for jv = 1:length(vars)
        jvar = (oo_.dr.obs_var == vars(jv));
        y = NaN(T,1); %Take care of T_eff instead of T for lags and NaN via mean with 'omitnan' option below
        y( (1-min(leadlags(jv),0)) : (T-max(leadlags(jv),0)), 1) = data( (1+max(leadlags(jv),0)) : (T+min(leadlags(jv),0)), jvar).^powers(jv);
        if jv==1
            m_data_tmp = y;
        else
            m_data_tmp = m_data_tmp.*y;
        end
    end
    % We replace NaN (due to leads and lags and missing values) with the corresponding mean
    if isoctave || matlab_ver_less_than('8.5')
        dataMoments(jm,1) = nanmean(m_data_tmp);
    else
        dataMoments(jm,1) = mean(m_data_tmp,'omitnan');
    end
    m_data_tmp(isnan(m_data_tmp)) = dataMoments(jm,1);
    m_data(:,jm) = m_data_tmp;
end


end %function end