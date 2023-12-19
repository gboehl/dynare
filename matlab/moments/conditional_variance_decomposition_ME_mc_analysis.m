function oo_ = ...
    conditional_variance_decomposition_ME_mc_analysis(NumberOfSimulations, type, dname, fname, Steps, exonames, exo, var_list, endo, mh_conf_sig, oo_,options_)
%oo_ = conditional_variance_decomposition_ME_mc_analysis(NumberOfSimulations, type, dname, fname, Steps, exonames, exo, var_list, endo, mh_conf_sig, oo_,options_)
% This function analyses the (posterior or prior) distribution of the
% endogenous variables' conditional variance decomposition with measurement error.
%
% INPUTS
%   NumberOfSimulations     [integer]           scalar, number of simulations.
%   type                    [string]            'prior' or 'posterior'
%   dname                   [string]            directory name where to save
%   fname                   [string]            name of the mod-file
%   Steps                   [integers]          horizons at which to conduct decomposition
%   exonames                [string]            (n_exo*char_length) character array with names of exogenous variables
%   exo                     [string]            name of current exogenous
%                                               variable
%   var_list                [string]            (n_endo*char_length) character array with name
%                                               of endogenous variables
%   endo                    [integer]           Current endogenous variable
%   mh_conf_sig             [double]            2 by 1 vector with upper
%                                               and lower bound of HPD intervals
%   oo_                     [structure]         Dynare structure where the results are saved.
%
% OUTPUTS
%   oo_          [structure]        Dynare structure where the results are saved.

% Copyright © 2017-2023 Dynare Team
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

if strcmpi(type,'posterior')
    TYPE = 'Posterior';
    PATH = [dname '/metropolis/'];
else
    TYPE = 'Prior';
    PATH = [dname '/prior/moments/'];
end

endogenous_variable_index = check_name(var_list, endo);
if isempty(endogenous_variable_index)
    disp([ type '_analysis:: Can''t find ' endo '!'])
    return
end

exogenous_variable_index = check_name(exonames,exo);
if isempty(exogenous_variable_index)
    if isequal(exo,'ME')
        exogenous_variable_index=length(exonames)+1;
    else
        disp([ type '_analysis:: ' exo ' is not a declared exogenous variable!'])
        return
    end
end

if isoctave && octave_ver_less_than('8.4') %Octave bug #60347
    [~,index_subset]=intersect_stable(var_list,options_.varobs);
else
    [~,index_subset]=intersect(var_list,options_.varobs,'stable');
end

matrix_pos=strmatch(endo, var_list(index_subset),'exact');
name_1 = endo;
name_2 = exo;
name = [ name_1 '.' name_2 ];

if isfield(oo_, [ TYPE 'TheoreticalMoments' ])
    temporary_structure = oo_.([TYPE 'TheoreticalMoments']);
    if isfield(temporary_structure,'dsge')
        temporary_structure = oo_.([TYPE 'TheoreticalMoments']).dsge;
        if isfield(temporary_structure,'ConditionalVarianceDecompositionME')
            temporary_structure = oo_.([TYPE 'TheoreticalMoments']).dsge.ConditionalVarianceDecompositionME.Mean;
            if isfield(temporary_structure,name)
                if sum(Steps-temporary_structure.(name)(1,:)) == 0
                    % Nothing (new) to do here...
                    return
                end
            end
        end
    end
end

ListOfFiles = dir([ PATH  fname '_' TYPE 'ConditionalVarianceDecompME*.mat']);
i1 = 1; tmp = zeros(NumberOfSimulations,length(Steps));
for file = 1:length(ListOfFiles)
    load([ PATH ListOfFiles(file).name ],'Conditional_decomposition_array_ME');
    % 4D-array (endovar,time,exovar,simul)
    i2 = i1 + size(Conditional_decomposition_array_ME,4) - 1;
    tmp(i1:i2,:) = transpose(dynare_squeeze(Conditional_decomposition_array_ME(matrix_pos,:,exogenous_variable_index,:)));
    i1 = i2+1;
end

p_mean = NaN(1,length(Steps));
p_median = NaN(1,length(Steps));
p_variance = NaN(1,length(Steps));
p_deciles = NaN(9,length(Steps));
if options_.estimation.moments_posterior_density.indicator
    p_density = NaN(2^9,2,length(Steps));
end
p_hpdinf = NaN(1,length(Steps));
p_hpdsup = NaN(1,length(Steps));
for i=1:length(Steps)
    if options_.estimation.moments_posterior_density.indicator
        [pp_mean, pp_median, pp_var, hpd_interval, pp_deciles, pp_density] = ...
            posterior_moments(tmp(:,i),1,mh_conf_sig);
        p_density(:,:,i) = pp_density;
    else
        [pp_mean, pp_median, pp_var, hpd_interval, pp_deciles] = ...
            posterior_moments(tmp(:,i),0,mh_conf_sig);
    end
    p_mean(i) = pp_mean;
    p_median(i) = pp_median;
    p_variance(i) = pp_var;
    p_deciles(:,i) = pp_deciles;
    p_hpdinf(i) = hpd_interval(1);
    p_hpdsup(i) = hpd_interval(2);
end

FirstField = sprintf('%sTheoreticalMoments', TYPE);

oo_.(FirstField).dsge.ConditionalVarianceDecompositionME.Steps = Steps;
oo_.(FirstField).dsge.ConditionalVarianceDecompositionME.Mean.(name_1).(name_2) = p_mean;
oo_.(FirstField).dsge.ConditionalVarianceDecompositionME.Median.(name_1).(name_2) = p_median;
oo_.(FirstField).dsge.ConditionalVarianceDecompositionME.Variance.(name_1).(name_2) = p_variance;
oo_.(FirstField).dsge.ConditionalVarianceDecompositionME.HPDinf.(name_1).(name_2) = p_hpdinf;
oo_.(FirstField).dsge.ConditionalVarianceDecompositionME.HPDsup.(name_1).(name_2) = p_hpdsup;
oo_.(FirstField).dsge.ConditionalVarianceDecompositionME.deciles.(name_1).(name_2)  = p_deciles;
if options_.estimation.moments_posterior_density.indicator
    oo_.(FirstField).dsge.ConditionalVarianceDecompositionME.density.(name_1).(name_2) = p_density;
end
