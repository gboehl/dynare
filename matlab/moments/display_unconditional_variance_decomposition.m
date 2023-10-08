function display_unconditional_variance_decomposition(M_,options_,oo_,ivar,stationary_vars,index_subset,ME_present)
% display_unconditional_variance_decomposition(M_,options_,oo_,ivar,stationary_vars,index_subset,ME_present)
% This function displays the unconditional variance decomposition 
%
% INPUTS
%   M_                  [struct]        Matlab's structure describing the Model (initialized by dynare, see @ref{M_}).          
%   options_            [struct]        Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%   oo_                 [struct]        structure describing the Model
%   i_var               [double]        Index of requested variables in declaration order
%   stationary_vars     [double]        index of stationary vars in requested output
%   index_subset        [integer]       index of observables in requested variables
%   ME_present          [boolean]       indicator whether measurement error is present for requested variables
%
% OUTPUTS
%   None

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

if M_.exo_nbr > 1
    skipline()
    if options_.order == 2
        title = 'APPROXIMATED VARIANCE DECOMPOSITION (in percent)';
    else
        title = 'VARIANCE DECOMPOSITION (in percent)';
    end
    title = add_filter_subtitle(title, options_);
    headers = M_.exo_names;
    headers = vertcat(' ', headers);
    labels=get_labels_transformed_vars(M_.endo_names,ivar(stationary_vars),options_,false);
    lh = cellofchararraymaxlength(labels)+2;
    dyntable(options_, title, headers, labels, 100*oo_.gamma_y{options_.ar+2}(stationary_vars,:), lh, 8, 2);
    if ME_present
        if isoctave && octave_ver_less_than('8.4') %Octave bug #60347
            [stationary_observables, pos_index_subset] = intersect_stable(index_subset, stationary_vars);
        else
            [stationary_observables, pos_index_subset] = intersect(index_subset, stationary_vars, 'stable');
        end
        headers_ME = vertcat(headers, 'ME');
        labels=get_labels_transformed_vars(M_.endo_names,ivar(stationary_observables),options_,false);
        dyntable(options_, [title,' WITH MEASUREMENT ERROR'], headers_ME, labels, ...
            oo_.variance_decomposition_ME(pos_index_subset,:), lh, 8, 2);
    end
    if options_.TeX
        headers = M_.exo_names_tex;
        headers = vertcat(' ', headers);
        labels=get_labels_transformed_vars(M_.endo_names_tex,ivar(stationary_vars),options_,true);
        lh = cellofchararraymaxlength(labels)+2;
        dyn_latex_table(M_, options_, title, 'th_var_decomp_uncond', headers, labels, 100*oo_.gamma_y{options_.ar+2}(stationary_vars,:), lh, 8, 2);
        if ME_present
            headers_ME = vertcat(headers, 'ME');
            labels=get_labels_transformed_vars(M_.endo_names_tex,ivar(stationary_observables),options_,true);
            dyn_latex_table(M_, options_, [title,' WITH MEASUREMENT ERROR'], ...
                'th_var_decomp_uncond_ME', headers_ME, labels, oo_.variance_decomposition_ME(pos_index_subset,:), lh, 8, 2);
        end
    end
end