function display_comparison_moments_irfs(M_, options_mom_, data_moments, model_moments)
% display_comparison_moments_irfs(M_, options_mom_, data_moments, model_moments)
% -------------------------------------------------------------------------
% Displays and saves to disk the comparison of the data moments/IRFs and the model moments/IRFs
% -------------------------------------------------------------------------
% INPUTS
% M_:             [structure]  model information
% options_mom_:   [structure]  method of moments options
% data_moments:   [vector]     data moments
% model_moments:  [vector]     model moments
% -------------------------------------------------------------------------
% OUTPUT
% No output, just displays and saves to disk the comparison of the data moments and the model moments
% -------------------------------------------------------------------------
% This function is called by
%  o mom.run
% -------------------------------------------------------------------------
% This function calls
% o dyn_latex_table
% o dyntable
% o cellofchararraymaxlength
% -------------------------------------------------------------------------

% Copyright © 2023 Dynare Team
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


if strcmp(options_mom_.mom.mom_method,'IRF_MATCHING')
    titl = upper('Comparison of matched data IRFs and model IRFs');
    headers = {'IRF','Data','Model'};
    idx = 1;
    for jj = 1:size(M_.matched_irfs,1)
        irf_varname = M_.matched_irfs{jj,1};
        irf_shockname = M_.matched_irfs{jj,2};
        % note that periods can span over multiple rows
        IRF_PERIODS = [];
        for kk = 1:size(M_.matched_irfs{jj,3},1)
            irf_periods = M_.matched_irfs{jj,3}{kk,1};
            IRF_PERIODS = [IRF_PERIODS; irf_periods(:)];
        end
        for hh = 1:length(IRF_PERIODS)
            labels{idx,1} = sprintf('%s %s (%u)',irf_varname,irf_shockname,IRF_PERIODS(hh));
            labels_TeX{idx,1} = sprintf('%s %s (%u)',M_.endo_names_tex{ismember(M_.endo_names,irf_varname)},M_.exo_names_tex{ismember(M_.exo_names,irf_shockname)},IRF_PERIODS(hh));
            idx = idx+1;
        end
    end
else
    titl = ['Comparison of matched data moments and model moments (',options_mom_.mom.mom_method,')'];
    headers = {'Moment','Data','Model'};
    for jm = 1:size(M_.matched_moments,1)
        lables_tmp = 'E[';
        lables_tmp_tex = 'E \left[ ';
        for jvar = 1:length(M_.matched_moments{jm,1})
            lables_tmp = [lables_tmp M_.endo_names{M_.matched_moments{jm,1}(jvar)}];
            lables_tmp_tex = [lables_tmp_tex, '{', M_.endo_names_tex{M_.matched_moments{jm,1}(jvar)}, '}'];
            if M_.matched_moments{jm,2}(jvar) ~= 0
                lables_tmp = [lables_tmp, '(', num2str(M_.matched_moments{jm,2}(jvar)), ')'];
                lables_tmp_tex = [lables_tmp_tex, '_{t', num2str(M_.matched_moments{jm,2}(jvar)), '}'];
            else
                lables_tmp_tex = [lables_tmp_tex, '_{t}'];
            end
            if M_.matched_moments{jm,3}(jvar) > 1
                lables_tmp = [lables_tmp, '^', num2str(M_.matched_moments{jm,3}(jvar))];
                lables_tmp_tex = [lables_tmp_tex, '^{', num2str(M_.matched_moments{jm,3}(jvar)) '}'];
            end
            if jvar == length(M_.matched_moments{jm,1})
                lables_tmp = [lables_tmp, ']'];
                lables_tmp_tex = [lables_tmp_tex, ' \right]'];
            else
                lables_tmp = [lables_tmp, '*'];
                lables_tmp_tex = [lables_tmp_tex, ' \times '];
            end
        end
        labels{jm,1} = lables_tmp;
        labels_TeX{jm,1} = lables_tmp_tex;
    end
end
data_mat = [data_moments model_moments];
dyntable(options_mom_, titl, headers, labels, data_mat, cellofchararraymaxlength(labels)+2, 10, 7);
if options_mom_.TeX
    dyn_latex_table(M_, options_mom_, titl, ['comparison_moments_', options_mom_.mom.mom_method], headers, labels_TeX, data_mat, cellofchararraymaxlength(labels)+2, 10, 7);
end