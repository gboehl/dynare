function smoother2histval(opts)
% This function takes values from oo_.SmoothedVariables (and possibly
% oo_.SmoothedShocks) and copies them into M_.histval.
%
% Optional fields in 'opts' structure:
%    infile:      An optional *_results MAT file created by Dynare.
%                 If present, oo_.Smoothed{Variables,Shocks} are read from
%                 there. Otherwise, they are read from the global workspace.
%    invars:      An optional char or cell array listing variables to read in
%                 oo_.SmoothedVariables. If absent, all the endogenous
%                 variables present in oo_.SmoothedVariables are used.
%    period:      An optional period number to use as the starting point
%                 for subsequent simulations. It should be between 1 and
%                 the number of observations that were used to produce the
%                 smoothed values. If absent, the last observation is used.
%    outfile:     An optional MAT file in which to save the histval structure.
%                 If absent, the output will be written in M_.endo_histval
%    outvars:     An optional char or cell array listing variables to be written in
%                 outfile or M_.endo_histval. This cell must be of same
%                 length than invars, and there is a mapping between the input
%                 variable at the i-th position in invars, and the output
%                 variable at the i-th position in outvars. If absent, then
%                 taken as equal to invars.
%
% The function also uses the value of option_.parameter_set

% Copyright (C) 2014-2021 Dynare Team
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

global M_ options_ oo_

if ~isfield(opts, 'infile')
    if ~isfield(oo_, 'SmoothedVariables')
        error('Could not find smoothed variables; did you set the "smoother" option?')
    end
    smoothedvars = oo_.SmoothedVariables;
    smoothedshocks = oo_.SmoothedShocks;
else
    S = load(opts.infile);
    if ~isfield(S, 'oo_') || ~isfield(S.oo_, 'SmoothedVariables')
        error('Could not find smoothed variables in file; is this a Dynare results file, and did you set the "smoother" option when producing it?')
    end
    smoothedvars = S.oo_.SmoothedVariables;
    smoothedshocks = S.oo_.SmoothedShocks;
end

% Hack to determine if oo_.SmoothedVariables was computed after a Metropolis
tmp = fieldnames(smoothedvars);
if isstruct(smoothedvars.(tmp{1}))
    post_metropolis = 1;
    if ~ isstruct(smoothedvars.(tmp{end}))
        % point and metropolis results are simultaneously present
        post_metropolis = 2;
    end

elseif isstruct(smoothedvars.(tmp{end}))
    % point and metropolis results are simultaneously present
    post_metropolis = 2;
else
    post_metropolis = 0;
end

if post_metropolis
    tmp = fieldnames(smoothedvars.Mean);
    if length(tmp)~=M_.endo_nbr
        warning(['You are using smoother2histval although smoothed values have not '...
            'been computed for all endogenous and auxiliary variables.'...
            'The value of these variables will be set to their steady state.'])
    end
    tmpexo = fieldnames(smoothedshocks.Mean);
else
    tmp = fieldnames(smoothedvars);
    tmpexo = fieldnames(smoothedshocks);
end

% If post-Metropolis, select the parameter set
if isempty(options_.parameter_set)
    if post_metropolis
        smoothedvars = smoothedvars.Mean;
        smoothedshocks = smoothedshocks.Mean;
    end
else
    switch options_.parameter_set
      case 'calibration'
        if post_metropolis == 1
            error('Option parameter_set=calibration is not consistent with computed smoothed values.')
        end
      case 'posterior_mode'
        if post_metropolis == 1
            error('Option parameter_set=posterior_mode is not consistent with computed smoothed values.')
        end
      case 'posterior_mean'
        if ~post_metropolis
            error('Option parameter_set=posterior_mean is not consistent with computed smoothed values.')
        end
        smoothedvars = smoothedvars.Mean;
        smoothedshocks = smoothedshocks.Mean;
      case 'posterior_median'
        if ~post_metropolis
            error('Option parameter_set=posterior_median is not consistent with computed smoothed values.')
        end
        smoothedvars = smoothedvars.Median;
        smoothedshocks = smoothedshocks.Median;
      otherwise
        error([ 'Option parameter_set=' options_.parameter_set ' unsupported.' ])
    end
end

% Determine number of periods
n = length(smoothedvars.(tmp{1}));

if n < M_.maximum_endo_lag
    error('Not enough observations to create initial conditions')
end

if isfield(opts, 'invars')
    invars = opts.invars;
    if ischar(invars)
        invars = cellstr(invars);
    end
else
    invars = [tmp; tmpexo];
end

if isfield(opts, 'period')
    period = opts.period;
    if period > n
        error('The period that you indicated is beyond the data sample')
    end
    if period < M_.maximum_endo_lag
        error('The period that you indicated is too small to construct initial conditions')
    end
else
    period = n;
end

if isfield(opts, 'outvars')
    outvars = opts.outvars;
    if ischar(outvars)
        outvars = cellstr(outvars);
    end
    if length(invars) ~= length(outvars)
        error('The number of input and output variables is not the same')
    end
else
    outvars = invars;
end

% Initialize outputs
if ~isfield(opts, 'outfile')
    % Output to M_.endo_histval
    M_.endo_histval = repmat(oo_.steady_state, 1, M_.maximum_lag);
else
    % Output to a file
    data = zeros(M_.maximum_endo_lag, length(invars));
    for i=1:length(outvars)
        j = strmatch(outvars{i}, M_.endo_names, 'exact');
        if ~isempty(j)
            data(:,i)=oo_.steady_state(j);
        end
    end
    o = dseries();
end

% Handle all endogenous variables to be copied
for i = 1:length(invars)
    if ~isempty(strmatch(invars{i}, M_.endo_names, 'exact')) 
        if oo_.Smoother.loglinear
            s = exp(smoothedvars.(invars{i}));
        else
            s = smoothedvars.(invars{i});
        end
    elseif ~isempty(strmatch(invars{i}, M_.exo_names, 'exact'))
        s = smoothedshocks.(invars{i});
    else
        error('smoother2histval: unknown input variable')
    end
    
    v = s((period-M_.maximum_lag+1):period);
    if ~isfield(opts, 'outfile')
        j_endo = strmatch(outvars{i}, M_.endo_names, 'exact');
        if ~isempty(j_endo) 
            M_.endo_histval(j_endo, :) = v;
        end
        j_exo = strmatch(outvars{i}, M_.exo_names, 'exact');
        if ~isempty(j_exo)
            M_.exo_histval(j_exo, :) = v;
        end
        if isempty(j_endo) && isempty(j_exo)
            error(['smoother2histval: output variable ' outvars{i} ' does not exist.'])            
        end
    else
        data(:, i) = v';
    end
end
if isfield(opts, 'outfile')
    o = dseries(data, '1Y', invars);
end

% $$$ % Handle auxiliary variables for lags (both on endogenous and exogenous)
% $$$ for i = 1:length(M_.aux_vars)
% $$$     if ~ ismember(M_.endo_names{M_.aux_vars(i).endo_index},invars)
% $$$         if M_.aux_vars(i).type ~= 1 && M_.aux_vars(i).type ~= 3
% $$$             continue
% $$$         end
% $$$         if M_.aux_vars(i).type == 1
% $$$             % Endogenous
% $$$             orig_var = M_.endo_names{M_.aux_vars(i).orig_index};
% $$$         else
% $$$             % Exogenous
% $$$             orig_var = M_.exo_names{M_.aux_vars(i).orig_index};
% $$$         end
% $$$         [m, k] = ismember(orig_var, outvars);
% $$$         if m
% $$$             if ~isempty(strmatch(invars{k}, M_.endo_names))
% $$$                 s = smoothedvars.(invars{k});
% $$$             else
% $$$                 s = smoothedshocks.(invars{k});
% $$$             end
% $$$             l = M_.aux_vars(i).orig_lead_lag;
% $$$             if period-M_.maximum_endo_lag+1+l < 1
% $$$                 error('The period that you indicated is too small to construct initial conditions')
% $$$             end
% $$$             j = M_.aux_vars(i).endo_index;
% $$$             v = s((period-M_.maximum_endo_lag+1+l):(period+l)); %+steady_state(j);
% $$$             if ~isfield(opts, 'outfile')
% $$$                 M_.endo_histval(j, :) = v;
% $$$             else
% $$$                 % When saving to a file, x(-2) is in the variable called "x_l2"
% $$$                 lead_lag = num2str(l);
% $$$                 lead_lag = regexprep(lead_lag, '-', 'l');
% $$$                 o.([ orig_var '_' lead_lag ]) = v;
% $$$             end
% $$$         end
% $$$     end
% $$$ end

% Finalize output
if isfield(opts, 'outfile')
    [dir, fname, ext] = fileparts(opts.outfile);
    if ~strcmp(ext,'.mat') && ~isempty(ext)
        error(['smoother2hisvtval: if outfile has an extension, it must ' ...
               'be .mat'])
    end
    o.save([dir fname]);
end

end
