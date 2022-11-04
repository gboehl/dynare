function [M_, options_] = setup(M_,options_, options_occbin_)
% function [M_, options_] = setup(M_, options_, options_occbin_)
% Sets up run of Occbin: creates shock matrix, sets options
%
% INPUT:
% - M_                  [structure]     Matlab's structure describing the model
% - options_            [structure]     Matlab's structure containing the options
% - options_occbin_     [structure]     Matlab's structure containing Occbin options
%
% OUTPUT:
% - M_                  [structure]     Matlab's structure describing the model
% - options_occbin_     [structure]     Matlab's structure containing Occbin options

% Copyright © 2021 Dynare Team
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

options_ = occbin.set_option(options_,options_occbin_,'simul.periods');
options_ = occbin.set_option(options_,options_occbin_,'simul.curb_retrench');
options_ = occbin.set_option(options_,options_occbin_,'simul.maxit');
options_ = occbin.set_option(options_,options_occbin_,'simul.check_ahead_periods');
options_ = occbin.set_option(options_,options_occbin_,'simul.debug');
options_ = occbin.set_option(options_,options_occbin_,'simul.periodic_solution');
options_ = occbin.set_option(options_,options_occbin_,'smoother.periods');
options_ = occbin.set_option(options_,options_occbin_,'smoother.curb_retrench');
options_ = occbin.set_option(options_,options_occbin_,'smoother.maxit');
options_ = occbin.set_option(options_,options_occbin_,'smoother.check_ahead_periods');
options_ = occbin.set_option(options_,options_occbin_,'smoother.debug');
options_ = occbin.set_option(options_,options_occbin_,'smoother.periodic_solution');
options_ = occbin.set_option(options_,options_occbin_,'smoother.inversion_filter');
options_ = occbin.set_option(options_,options_occbin_,'likelihood.periods');
options_ = occbin.set_option(options_,options_occbin_,'likelihood.curb_retrench');
options_ = occbin.set_option(options_,options_occbin_,'likelihood.maxit');
options_ = occbin.set_option(options_,options_occbin_,'likelihood.check_ahead_periods');
options_ = occbin.set_option(options_,options_occbin_,'likelihood.periodic_solution');
options_ = occbin.set_option(options_,options_occbin_,'likelihood.max_number_of_iterations');
options_ = occbin.set_option(options_,options_occbin_,'filter.use_relaxation');
options_ = occbin.set_option(options_,options_occbin_,'likelihood.inversion_filter');

if isfield(M_,'surprise_shocks') && ~isempty(M_.surprise_shocks)
    temp=zeros(max(cat(2,M_.surprise_shocks.periods)),M_.exo_nbr);
    for ii = 1:length(M_.surprise_shocks)
        ivar = M_.surprise_shocks(ii).exo_id;
        temp(M_.surprise_shocks(ii).periods,ivar) = M_.surprise_shocks(ii).value;
    end
    shock_index=~all(temp==0,1);
    options_.occbin.simul.SHOCKS=temp(:,shock_index);
    options_.occbin.simul.exo_pos=find(shock_index);
end
