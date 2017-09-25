function [dbase, info] = checkdatabaseforinversion(dbase, DynareModel)

% Check that dbase contains all the endogenous variables of the model, and
% reorder the endogenous variables as declared in the mod file. If Dynare
% adds auxiliary variables, for lags greater than 1 on endogebnous variables
% or lags on the exogenous variables.

% Copyright (C) 2017 Dynare Team
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

listoflaggedexogenousvariables = {};

info = struct;

k = 0;
for i = DynareModel.orig_endo_nbr+1:DynareModel.endo_nbr
    k = k+1;
    if DynareModel.aux_vars(k).type==1
        if ismember(deblank(DynareModel.endo_names(DynareModel.aux_vars(k).orig_index,:)), dbase.name)
            dbase{deblank(DynareModel.endo_names(DynareModel.aux_vars(k).endo_index, :))} = dbase{deblank(DynareModel.endo_names(DynareModel.aux_vars(k).orig_index, :))}.lag(abs(DynareModel.aux_vars(k).orig_lead_lag));
        else
            error('%s not available in dbase!', deblank(DynareModel.endo_names(DynareModel.aux_vars(k).orig_index, :)));
        end
    elseif DynareModel.aux_vars(k).type==3 && DynareModel.aux_vars(k).orig_lead_lag==0
        dbase{deblank(DynareModel.endo_names(DynareModel.aux_vars(k).endo_index,:))} = dbase{deblank(DynareModel.exo_names(DynareModel.aux_vars(k).orig_index, :))}.lag(abs(DynareModel.aux_vars(k).orig_lead_lag));
        listoflaggedexogenousvariables = vertcat(listoflaggedexogenousvariables, deblank(DynareModel.exo_names(DynareModel.aux_vars(k).orig_index, :)));
    else
        warning('Please contact Dynare Team!')
    end
end

info.endonames = cellstr(DynareModel.endo_names);
info.exonames = cellstr(DynareModel.exo_names);
info.computeresiduals = false;

% Check that all the endogenous variables are defined in dbase.
missingendogenousvariables = setdiff(info.endonames, dbase.name);
if ~isempty(missingendogenousvariables)
    disp('Some endognous variables are missing:')
    missingendogenousvariables
    error()
end

% Check if all the exogenous variables are in dbase.
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