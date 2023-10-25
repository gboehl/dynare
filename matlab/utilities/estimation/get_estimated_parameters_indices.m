function ipnames = get_estimated_parameters_indices(params, pnames, eqname, M_)
% ipnames = get_estimated_parameters_indices(params, pnames, eqname, M_)

% Copyright © 2021-2023 Dynare Team
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

list_of_parameters = fieldnames(params);

% Check that the estimated parameters are used in the PAC equation.
parameters_not_in_equation = setdiff(list_of_parameters, pnames);
if ~isempty(parameters_not_in_equation)
    skipline()
    if length(parameters_not_in_equation)>1
        list = sprintf('  %s\n', parameters_not_in_equation{:});
        remk = sprintf('  The following parameters:\n\n%s\n  do not appear in equation %s.', list, eqname);
    else
        remk = sprintf('  Parameter %s does not appear in equation %s.', parameters_not_in_equation{1}, eqname);
    end
    disp(remk)
    skipline()
    error('The estimated parameters must be used in equation %s.', eqname)
end

ipnames = zeros(size(list_of_parameters));
for i=1:length(ipnames)
    ipnames(i) = find(strcmp(list_of_parameters{i}, M_.param_names));
end