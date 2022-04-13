function mh_autocorrelation_function(options_,M_,estim_params_,type,blck,name1,name2)
% This function plots the autocorrelation of the sampled draws in the
% posterior distribution.
%
%
% INPUTS
%
%   options_        [structure]    Dynare structure.
%   M_              [structure]    Dynare structure (related to model definition).
%   estim_params_   [structure]    Dynare structure (related to estimation).
%   type            [string]       'DeepParameter', 'MeasurementError' (for measurement equation error) or 'StructuralShock' (for structural shock).
%   blck            [integer]      Number of the mh chain.
%   name1           [string]       Object name.
%   name2           [string]       Object name.
%
% OUTPUTS
%   None
%
% SPECIAL REQUIREMENTS

% Copyright © 2003-2017 Dynare Team
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

% Cet the column index:
if nargin<7
    column = name2index(options_, M_, estim_params_, type, name1);
else
    column = name2index(options_, M_, estim_params_, type, name1, name2);
end

if isempty(column)
    return
end

% Get informations about the posterior draws:
MetropolisFolder = CheckPath('metropolis',M_.dname);
load_last_mh_history_file(MetropolisFolder, M_.fname);

FirstMhFile = record.KeepedDraws.FirstMhFile;
FirstLine = record.KeepedDraws.FirstLine; ifil = FirstLine;
TotalNumberOfMhFiles = sum(record.MhDraws(:,2));
TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
clear record;

% Get all the posterior draws:
PosteriorDraws = GetAllPosteriorDraws(column, FirstMhFile, FirstLine, TotalNumberOfMhFiles, NumberOfDraws, blck);

% Compute the autocorrelation function:
[autocov,autocor] = sample_autocovariance(PosteriorDraws,options_.mh_autocorrelation_function_size);

% Plot the posterior draws:

if strcmpi(type,'DeepParameter')
    TYPE = 'parameter ';
elseif strcmpi(type,'StructuralShock')
    if nargin<7
        TYPE = 'the standard deviation of structural shock ';
    else
        TYPE = 'the correlation between structural shocks ';
    end
elseif strcmpi(type,'MeasurementError')
    if nargin<7
        TYPE = 'the standard deviation of measurement error ';
    else
        TYPE = 'the correlation between measurement errors ';
    end
end

if nargin<7
    FigureName = ['Autocorrelogram for ' TYPE name1];
else
    FigureName = ['Autocorrelogram for ' TYPE name1 ' and ' name2];
end

if options_.mh_nblck>1
    FigureName = [ FigureName , ' (block number' int2str(blck)  ').'];
end

hh=dyn_figure(options_.nodisplay,'Name',FigureName);

bar(0:options_.mh_autocorrelation_function_size,autocor,'k');
axis tight
% create subdirectory <dname>/graphs if it doesn't exist
if ~exist(M_.dname, 'dir')
    mkdir('.',M_.dname);
end
if ~exist([M_.dname filesep 'graphs'])
    mkdir(M_.dname,'graphs');
end

plot_name=get_the_name(column,0,M_,estim_params_,options_);
dyn_saveas(hh,[M_.dname, filesep, 'graphs', filesep, 'MH_Autocorrelation_' plot_name],options_.nodisplay,options_.graph_format)
