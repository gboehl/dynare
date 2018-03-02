function [dbase, info] = checkdatabase(dbase, DynareModel, inversionflag)

% Check that dbase contains all the endogenous variables of the model, and
% reorder the endogenous variables as declared in the mod file. If Dynare
% adds auxiliary variables, for lags greater than 1 on endogenous variables,
% endogenous variables in difference (which may be lagged), or lags on the
% exogenous variables, then thee routine complete the database.

% Copyright (C) 2018 Dynare Team
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

if DynareModel.maximum_endo_lead
    error('The model (%s) is assumed to be backward!', DynareModel.fname)
end

if nargin<3
    inversionflag = false;
end

listoflaggedexogenousvariables = {};

info = struct;

k = 0;
for i = DynareModel.orig_endo_nbr+1:DynareModel.endo_nbr
    k = k+1;
    if DynareModel.aux_vars(k).type==1
        if ismember(DynareModel.endo_names{DynareModel.aux_vars(k).orig_index}, dbase.name)
            dbase{DynareModel.endo_names{DynareModel.aux_vars(k).endo_index}} = dbase{DynareModel.endo_names{DynareModel.aux_vars(k).orig_index}}.lag(abs(DynareModel.aux_vars(k).orig_lead_lag));
        else
            error('%s not available in dbase!', DynareModel.endo_names{DynareModel.aux_vars(k).orig_index});
        end
    elseif DynareModel.aux_vars(k).type==3
        dbase{DynareModel.endo_names{DynareModel.aux_vars(k).endo_index}} = dbase{DynareModel.exo_names{DynareModel.aux_vars(k).orig_index}}.lag(abs(DynareModel.aux_vars(k).orig_lead_lag));
        listoflaggedexogenousvariables = vertcat(listoflaggedexogenousvariables, DynareModel.exo_names{DynareModel.aux_vars(k).orig_index});
    elseif DynareModel.aux_vars(k).type==8
        dbase{DynareModel.endo_names{DynareModel.aux_vars(k).endo_index}} = dbase{DynareModel.endo_names{DynareModel.aux_vars(k).orig_index}}.diff.lag(abs(DynareModel.aux_vars(k).orig_lead_lag));
    else
        warning('Please contact Dynare Team!')
    end
end

info.endonames = DynareModel.endo_names;
info.exonames = DynareModel.exo_names;
info.computeresiduals = false;

% Check that all the endogenous variables are defined in dbase.
missingendogenousvariables = setdiff(info.endonames, dbase.name);
if ~isempty(missingendogenousvariables)
    disp('Some endognous variables are missing:')
    missingendogenousvariables
    error()
end

if inversionflag
    % If some exogenous variables are missing, check that they can be
    % interpreted as residuals. 
    missingexogenousvariables = setdiff(info.exonames, dbase.name);
    if ~isempty(missingexogenousvariables)
        disp(sprintf('%s exogenous variables are missing in the database...', num2str(length(missingexogenousvariables))))
        listofmissinglaggedexognousvariables = intersect(listoflaggedexogenousvariables, missingexogenousvariables);
        if isempty(listofmissinglaggedexognousvariables)
            info.residuals = missingexogenousvariables;
            info.computeresiduals = true;
            disp('These variables can be calibrated by calling calibrateresiduals routine.')
        else
            info.residuals = setdiff(missingexogenousvariables, listofmissinglaggedexognousvariables);
            disp('The following exogenous variables:')
            listofmissinglaggedexognousvariables
            disp('are not residuals, and cannot be calibrated by calling calibrateresiduals.')
        end
    else
        disp('All the endogenous and exogenous variables are calibrated!')
    end
else
    % Check that all the exogenous variables are defined in dbase
    missingexogenousvariables = setdiff(info.exonames, dbase.name);
    if ~isempty(missingexogenousvariables)
        disp('Some exognous variables are missing:')
        missingexogenousvariables
        error()
    end
    info.residuals = [];
end