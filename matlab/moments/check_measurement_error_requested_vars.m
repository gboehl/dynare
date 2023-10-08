function [ME_present,observable_pos_requested_vars,index_subset,index_observables]=check_measurement_error_requested_vars(M_,options_,ivar)
% [ME_present,observable_pos_requested_vars,index_subset,index_observables]=check_measurement_error_requested_vars(M_,options_,ivar)
% This function checks for the presence of measurement error and outputs
% the indices of the affected variables
%
% INPUTS
%   M_                  [struct]        Matlab's structure describing the Model
%   options_            [struct]        Matlab's structure describing the options
%   i_var :             [double]        Index of requested variables in declaration order
%
% OUTPUTS
%   ME_present                      [boolean]       indicator whether measurement error is present for requested variables
%   observable_pos_requested_vars   [integer]       index of observables in list of endogenous variables
%   index_subset                    [integer]       index of observables in ivar
%   index_observables               [integer]       index of requested i_var in observables 

% Copyright © 2023 Dynare Team
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


index_subset=[];
index_observables=[];
observable_pos_requested_vars=[];

ME_present=false;
if ~all(diag(M_.H)==0)
    if isoctave && octave_ver_less_than('8.4') %Octave bug #60347
        [observable_pos_requested_vars,index_subset,index_observables]=intersect_stable(ivar,options_.varobs_id);
    else
        [observable_pos_requested_vars,index_subset,index_observables]=intersect(ivar,options_.varobs_id,'stable');
    end
    if ~isempty(observable_pos_requested_vars)
        ME_present=true;
    end
end
