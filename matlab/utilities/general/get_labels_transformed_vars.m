function labels=get_labels_transformed_vars(endo_names,var_indices,options_,TeX)
% function labels=get_labels_transformed_vars(endo_names,var_indices,options_,TeX)
% This function provides the variable labels for table outputs in case of
% applied transformations like logs
%
% INPUTS
%   endo_names          [cell]        cell array of variable names
%   var_indices         [double]      vector of variable indices
%   options_            [structure]   Dynare structure containing the options
%   TeX                 [boolean]     indicator for TeX-output
% OUTPUTS
%   labels              [cell]        cell array of variable labels
%
% Copyright © 2020 Dynare Team
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

if nargin<4
    TeX=0;
end
if options_.loglinear
    labels=[];
    for var_iter=1:length(var_indices)
        if TeX
            labels{var_iter,1}=['\log(',endo_names{var_indices(var_iter)},')'];
        else
            labels{var_iter,1}=['log(',endo_names{var_indices(var_iter)},')'];
        end
    end
else
    labels = endo_names(var_indices);
end
end