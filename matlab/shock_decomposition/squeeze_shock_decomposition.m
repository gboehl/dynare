function oo_ = squeeze_shock_decomposition(M_,oo_,options_,var_list_)
%function oo_ = squeeze_shock_decomposition(M_,oo_,options_,var_list_)
% INPUTS
%    M_:          [structure]            Definition of the model
%    oo_:         [structure]            Storage of results
%    options_:    [structure]            Options
%    var_list_:   [cell of char arrays]  List of variables
%
% OUTPUTS
%    oo_:         [structure]            Storage of results
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2019 Dynare Team
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

if ~options_.shock_decomp.with_epilogue
    endo_names = M_.endo_names;
else
    endo_names = [M_.endo_names; M_.epilogue_names];

end
if isfield(oo_,'plot_shock_decomposition_info') && isfield(oo_.plot_shock_decomposition_info','i_var')
    my_vars = oo_.plot_shock_decomposition_info.i_var;
else
    my_vars=[];
end
if nargin>3
    my_vars = [varlist_indices(var_list_,endo_names); my_vars];
end
sd_vlist = endo_names(my_vars,:);

if isfield(options_.plot_shock_decomp,'q2a') && isstruct(options_.plot_shock_decomp.q2a)

    avname={options_.plot_shock_decomp.q2a.qname};
    sda = options_.plot_shock_decomp.q2a(ismember(avname,sd_vlist));
    for k=1:length(sda)
        if isstruct(sda(k).aux)
            sd_vlist = [sd_vlist; {sda(k).aux.y}];
        end
    end
end

if isempty(sd_vlist)
    disp('Nothing has been squeezed: there is no list of variables for it!')
    return
end
i_var = varlist_indices(sd_vlist,endo_names);

oo_.plot_shock_decomposition_info.i_var = i_var;
oo_.shock_decomposition_info.i_var = i_var;
if isfield (oo_,'shock_decomposition')
    oo_.shock_decomposition = oo_.shock_decomposition(i_var,:,:);
end
if isfield (oo_,'realtime_conditional_shock_decomposition')
    oo_.realtime_conditional_shock_decomposition = ...
        my_squeeze(oo_.realtime_conditional_shock_decomposition, i_var);
end
if isfield (oo_,'realtime_forecast_shock_decomposition')
    oo_.realtime_forecast_shock_decomposition = ...
        my_squeeze(oo_.realtime_forecast_shock_decomposition, i_var);
end
if isfield (oo_,'realtime_shock_decomposition')
    oo_.realtime_shock_decomposition = ...
        my_squeeze(oo_.realtime_shock_decomposition, i_var);
end
if isfield (oo_,'conditional_shock_decomposition')
    oo_.conditional_shock_decomposition = ...
        my_squeeze(oo_.conditional_shock_decomposition, i_var);
end
if isfield (oo_,'initval_decomposition')
    oo_.initval_decomposition = oo_.initval_decomposition(i_var,:,:);
end

end

function shock_decomposition = my_squeeze(shock_decomposition, i_var)
fnam = fieldnames(shock_decomposition);
for k=1:length(fnam)
    shock_decomposition.(fnam{k}) =  shock_decomposition.(fnam{k})(i_var,:,:);
end

end
