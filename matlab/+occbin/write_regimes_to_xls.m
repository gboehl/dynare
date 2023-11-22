function write_regimes_to_xls(occbin_struct,M_,options_)
% function write_regimes_to_xls(occbin_struct,M_,options_)
% writes regime results to Excel-file
%
% INPUTS
% - occbin_struct   [struct]    occbin structure containing information on the regimes
% - M_              [struct]    Matlab's structure describing the model
% - options_        [struct]    Matlab's structure describing the current options

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

OutputDirectoryName = CheckPath('Output',M_.dname);

if strcmpi(options_.occbin.write_regimes.type,'simul') || strcmpi(options_.occbin.write_regimes.type,'smoother')
    if isfield(occbin_struct,options_.occbin.write_regimes.type) && isfield(occbin_struct.(options_.occbin.write_regimes.type),'regime_history')
        regime_history=occbin_struct.(lower(options_.occbin.write_regimes.type)).regime_history;
    else
        error('write_regimes_to_xls: the required field does not exist');    
    end
else
    error('write_regimes_to_xls: output type can only be simul or smoother.')
end

if isempty(options_.occbin.write_regimes.periods)
    T=1:length(regime_history);
else
    T=options_.occbin.write_regimes.periods;
end

xls_filename = options_.occbin.write_regimes.filename;

if isfield(regime_history,'regime')
    Header = {'time', 'regime_sequence', 'starting_period_of_regime'};
    for tp=1:length(T)
        xlsmat{tp,1}=T(tp);
        xlsmat{tp,2}=int2str(regime_history(tp).regime);
        xlsmat{tp,3}=int2str(regime_history(tp).regimestart);
    end
else
    Header = {'time', 'regime_sequence_1', 'starting_period_of_regime_1', 'regime_sequence_2', 'starting_period_of_regime_2'};
    for tp=1:length(T)
        xlsmat{tp,1}=T(tp);
        xlsmat{tp,2}=int2str(regime_history(tp).regime1);
        xlsmat{tp,3}=int2str(regime_history(tp).regimestart1);
        xlsmat{tp,4}=int2str(regime_history(tp).regime2);
        xlsmat{tp,5}=int2str(regime_history(tp).regimestart2);
    end
end

filename=[OutputDirectoryName filesep xls_filename '.xls'];

if isfile(filename)
    delete(filename)
end

if isoctave
    % “writetable” and “array2table” don’t exist under Octave
    if isoctave && ~user_has_octave_forge_package('io')
        error('The io package is required to write XLS files from Octave')
    end
    xlswrite(filename, vertcat(Header, xlsmat));
else
    writetable(array2table(xlsmat,'VariableNames',Header), filename, 'Sheet', 'Regimes');
end
