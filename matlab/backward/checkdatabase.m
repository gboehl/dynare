function [dbase, info] = checkdatabase(dbase, M_, inversionflag, simulationflag)
% [dbase, info] = checkdatabase(dbase, M_, inversionflag, simulationflag)
% Check that dbase contains all the endogenous variables of the model, and
% reorder the endogenous variables as declared in the mod file. If Dynare
% adds auxiliary variables, for lags greater than 1 on endogenous variables,
% endogenous variables in difference (which may be lagged), or lags on the
% exogenous variables, then thee routine complete the database.

% Copyright © 2018-2023 Dynare Team
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

% if M_.maximum_endo_lead
%    error('The model (%s) is assumed to be backward!', M_.fname)
% end

if nargin<3
    inversionflag = false;
end

if exist(sprintf('+%s/dynamic_set_auxiliary_series.m', M_.fname), 'file')
    dbase = feval(sprintf('%s.dynamic_set_auxiliary_series', M_.fname), dbase, M_.params);
end

listoflaggedexogenousvariables = {};
if ~isempty(M_.aux_vars)
    listoflaggedexogenousvariables = M_.exo_names([M_.aux_vars(find([M_.aux_vars.type]==3)).orig_index]);
end

listoflaggedendogenousvariables = {};
laggedendogenousvariablesidx = find(M_.lead_lag_incidence(1,1:M_.orig_endo_nbr));
if ~isempty(laggedendogenousvariablesidx)
    listoflaggedendogenousvariables = M_.endo_names(laggedendogenousvariablesidx);
end
if ~isempty(M_.aux_vars)
    laggedendogenousvariablesidx = find([M_.aux_vars.type]==1);
    if ~isempty(laggedendogenousvariablesidx)
        listoflaggedendogenousvariables = union(listoflaggedendogenousvariables, M_.endo_names([M_.aux_vars(laggedendogenousvariablesidx).orig_index]));
    end
    laggedendogenousvariablesidx = find([M_.aux_vars.type]==8);
    if ~isempty(laggedendogenousvariablesidx)
        listoflaggedendogenousvariables = union(listoflaggedendogenousvariables, M_.endo_names([M_.aux_vars(laggedendogenousvariablesidx).orig_index]));
    end
end

info = struct;
info.endonames = M_.endo_names;
info.exonames = M_.exo_names;
info.computeresiduals = false;

% Check that all the endogenous variables are defined in dbase.
missingendogenousvariables = setdiff(info.endonames, dbase.name);
if ~isempty(missingendogenousvariables)
    missinglaggedendogenousvariables = intersect(missingendogenousvariables, listoflaggedendogenousvariables);
    if ~isempty(missinglaggedendogenousvariables)
        error('Endogenous variable %s is missing in the database!', missinglaggedendogenousvariables{:})
    end
end

if inversionflag
    % If some exogenous variables are missing, check that they can be interpreted as residuals.
    missingexogenousvariables = setdiff(info.exonames, dbase.name);
    if ~isempty(missingexogenousvariables)
        dprintf('%s exogenous variables are missing in the database...', num2str(length(missingexogenousvariables)))
        listofmissinglaggedexognousvariables = intersect(listoflaggedexogenousvariables, missingexogenousvariables);
        if isempty(listofmissinglaggedexognousvariables)
            info.residuals = missingexogenousvariables;
            info.computeresiduals = true;
            disp('These variables can be calibrated by calling calibrateresiduals routine.')
        else
            info.residuals = setdiff(missingexogenousvariables, listofmissinglaggedexogenousvariables);
            disp('The following exogenous variables:')
            listofmissinglaggedexognousvariables
            disp('are not residuals, and cannot be calibrated by calling calibrateresiduals.')
        end
    else
        disp('All the endogenous and exogenous variables are calibrated!')
    end
elseif simulationflag
    % Check that all the exogenous variables are defined in dbase
    missingexogenousvariables = setdiff(info.exonames, dbase.name);
    listofmissingexovarforinit = intersect(missingexogenousvariables, listoflaggedexogenousvariables);
    if ~isempty(listofmissingexovarforinit)
        error('Exogenous variable %s is missing!', listofmissingexovarforinit{:})
    end
    info.residuals = [];
end
