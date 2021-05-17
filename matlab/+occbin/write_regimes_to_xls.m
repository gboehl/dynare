function write_regimes_to_xls(regime_history,M_,options_)
% function write_regimes_to_xls(regime_history,M_,options_)
% writes regime results to Excel-file
%
% INPUTS
% - regime_history  [struct]    information on the regimes
% - M_              [struct]    Matlab's structure describing the model
% - options_        [struct]    Matlab's structure describing the current options

% Copyright (C) 2021 Dynare Team
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

if isempty(options_.occbin.write_regimes.periods)
    T=1:length(regime_history);
else
    T=options_.occbin.write_regimes.periods;
end

xls_filename = options_.occbin.write_regimes.filename;

if isfield(regime_history,'regime')
    Header = {'time', 'regime sequence', 'starting period of regime'};
    for tp=1:length(T)
        xlsmat{tp,1}=T(tp);
        xlsmat{tp,2}=int2str(regime_history(tp).regime);
        xlsmat{tp,3}=int2str(regime_history(tp).regimestart);
    end
else
    Header = {'time', 'regime sequence 1', 'starting period of regime 1', 'regime sequence 2', 'starting period of regime 2'};
    for tp=1:length(T)
        xlsmat{tp,1}=T(tp);
        xlsmat{tp,2}=int2str(regime_history(tp).regime1);
        xlsmat{tp,3}=int2str(regime_history(tp).regimestart1);
        xlsmat{tp,4}=int2str(regime_history(tp).regime2);
        xlsmat{tp,5}=int2str(regime_history(tp).regimestart2);
    end
end
filename=[OutputDirectoryName filesep xls_filename '.xls'];
if matlab_ver_less_than('9.3')
    if exist(filename,'file')
        delete(filename)
    end    
else
    if isfile(filename)
        delete(filename)
    end
end
writetable(array2table(xlsmat,'VariableNames',Header), filename, 'Sheet', 'Regimes'); 
