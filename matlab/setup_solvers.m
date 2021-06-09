function Model = setup_solvers(Model)

% Setup solve_algo={12,14} by identifying equations with a log on the left hand side.
%
% INPUTS
% - Model     [struct]      Model description, aka M_.
% - Options   [struct]      Dynare's options, aka options_.
%
% OUTPUTS
% - Model     [struct]      Updated model description.

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

cannot_use_solve_algo_12_14 = false;

try
    json = loadjson_(sprintf('%s/model/json/modfile.json', Model.fname));
catch
    cannot_use_solve_algo_12_14 = true;
    message = 'Algorithms solve_algo={12,14} require json output of the model (use json=compute option)';
end

if ~cannot_use_solve_algo_12_14

    lhs = cell(length(json.model),1);
    isauxdiffloggedrhs = false(length(json.model), 1);

    for i = 1:length(json.model)
        if length(json.model)>1
            lhs{i} = json.model{i}.lhs;
        else
            lhs{i} = json.model.lhs;
        end
        if isempty(regexp(lhs{i}, '^\w+$|^log\(\w+\)$'))
            cannot_use_solve_algo_12_14 = true;
            message = sprintf('With solve_algo={12,14}, each equation must have on the left hand side a single variable or logged variable (equation %d does not satisfy this condition).', i);
            break
        end
        if length(json.model)>1
            rhs = json.model{i}.rhs;
        else
            rhs = json.model.lhs;
        end
        if i>Model.orig_endo_nbr &&  ~isempty(regexp(lhs{i}, '\<AUX_DIFF_(\d*)\>', 'once')) && ismember(lhs{i}, lhs(1:i-1)) && ...
                ~isempty(regexp(rhs, 'log\(\w*\)-log\(\w*\(-1\)\)', 'once'))
            isauxdiffloggedrhs(i) = true;
        end
    end

end

if ~cannot_use_solve_algo_12_14

    islog = @(x) ~isempty(regexp(x, 'log\(\w*\)', 'once'));
    
    lhs0 = lhs;
    for i=1:length(json.model)
        if islog(lhs{i})
            lhs0{i} = strrep(strrep(lhs{i}, 'log(', ''), ')', '');
        end
    end
    
    if ~isequal(length(unique(lhs0(1:Model.orig_endo_nbr))), length(lhs0(1:Model.orig_endo_nbr)))
        cannot_use_solve_algo_12_14 = true;
        message = sprintf('With solve_algo={12,14}, each equation must determine a different endogenous variable.')
    end

end

if cannot_use_solve_algo_12_14
    Model.isloggedlhs = {};
    Model.lhs = {};
    Model.isauxdiffloggedrhs = [];
    Model.possible_to_use_solve_algo_12_14 = false;
    Model.message_solve_algo_12_14 = message;
else
    Model.isloggedlhs = cellfun(islog, lhs);
    Model.lhs = lhs;
    Model.isauxdiffloggedrhs = isauxdiffloggedrhs;
    Model.possible_to_use_solve_algo_12_14 = true;
    Model.message_solve_algo_12_14 = '';
end