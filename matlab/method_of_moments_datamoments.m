function [dataMoments, m_data] = method_of_moments_datamoments(data, DynareResults, MatchedMoments, OptionsMoM)
% [dataMoments, m_data] = method_of_moments_datamoments(data, DynareResults, MatchedMoments, OptionsMoM)
% This function computes the user-selected empirical moments from data
% =========================================================================
% INPUTS
%  o data                    [T x varobs_nbr]  data set
%  o DynareResults:          [structure]       storage for results (oo_)
%  o MatchedMoments:         [structure]       information about selected moments to match in estimation (matched_moments_)
%  o OptionsMoM:             [structure]       information about all settings (specified by the user, preprocessor, and taken from global options_)
% -------------------------------------------------------------------------
% OUTPUTS
%  o dataMoments             [numMom x 1]       mean of selected empirical moments
%  o m_data                  [T x numMom]       selected empirical moments at each point in time
% -------------------------------------------------------------------------
% This function is called by
%  o method_of_moments.m
%  o method_of_moments_SMM.m
% =========================================================================
% Copyright (C) 2020 Dynare Team
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
% -------------------------------------------------------------------------
% Author(s): 
% o Willi Mutschler (willi@mutschler.eu)
% o Johannes Pfeifer (jpfeifer@uni-koeln.de)
% =========================================================================

% Initialization
T = size(data,1); % Number of observations (T) and number of observables (ny)
mom_nbr = OptionsMoM.mom_nbr;
dataMoments = nan(mom_nbr,1);
m_data = nan(T,mom_nbr);
% Product moment for each time period, i.e. each row t contains yt1(l1)^p1*yt2(l2)^p2*...
% note that here we already are able to treat leads and lags and any power product moments
for jm = 1:mom_nbr
    vars     = DynareResults.dr.inv_order_var(MatchedMoments{jm,1})';
    leadlags = MatchedMoments{jm,2}; % note that lags are negative numbers and leads are positive numbers
    powers   = MatchedMoments{jm,3};
    for jv = 1:length(vars)
        jvar = DynareResults.dr.obs_var == vars(jv);
        y = nan(T,1);
        y( (1-min(leadlags(jv),0)) : (T-max(leadlags(jv),0)) , 1) = data( (1+max(leadlags(jv),0)) : (T+min(leadlags(jv),0)) , jvar).^powers(jv);
        if jv==1
            m_data_tmp = y;
        else
            m_data_tmp = m_data_tmp.*y;
        end
    end
    dataMoments(jm,1) = sum(m_data_tmp,'omitnan')/(T-sum(abs(leadlags)));
    % We replace nan (due to leads and lags) with the corresponding mean
    % @wmutschl: this should also work for missing values, right?
    m_data_tmp(isnan(m_data_tmp)) = dataMoments(jm,1);
    m_data(:,jm) = m_data_tmp;
end


end %function end



