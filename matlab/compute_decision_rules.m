function [dr,info,M_,oo_] =compute_decision_rules(M_,options_,oo_)
% function [dr,info,M_,oo_] =compute_decision_rules(M_,options_,oo_)
% INPUTS
% - M_            [structure]     Matlab's structure describing the model (M_).
% - options_      [structure]     Matlab's structure describing the current options (options_).
% - oo_           [structure]     Matlab's structure containing the results (oo_).
%
% OUTPUTS
% - dr            [structure]     Reduced form model.
% - info          [integer]       scalar or vector, error code.
% - M_            [structure]     Matlab's structure describing the model (M_).
% - oo_           [structure]     Matlab's structure containing the results (oo_).

% Copyright (C) 2020 Dynare Team
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

if options_.discretionary_policy
    [dr,info,M_,oo_] = discretionary_policy_1(options_.instruments,M_,options_,oo_);
else
    [dr,info,M_,oo_] = resol(0,M_,options_,oo_);
end
