function check_varobs_are_endo_and_declared_once(varobs,endo_names)
% function check_varobs_are_endo_and_declared_once(varobs,endo_names)
% -------------------------------------------------------------------------
% Check that each declared observed variable:
% - is also an endogenous variable
% - is declared only once
% -------------------------------------------------------------------------
% INPUTS
%  o varobs:                 [cell] list of observed variables
%  o endo_names:             [cell] list of endogenous variables
% -------------------------------------------------------------------------
% OUTPUTS
%  none, display an error message something is wrong with VAROBS
% -------------------------------------------------------------------------
% This function is called by
%  o dynare_estimation_init.m
% -------------------------------------------------------------------------

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

number_of_observed_variables = length(varobs);
for i = 1:number_of_observed_variables    
    if ~any(strcmp(varobs{i},endo_names))
        error(['VAROBS: unknown variable (' varobs{i} ')!'])
    end
end

% Check that a variable is not declared as observed more than once.
if length(unique(varobs))<length(varobs)
    for i = 1:number_of_observed_variables
        if sum(strcmp(varobs{i},varobs)) > 1        
            error(['VAROBS: a variable cannot be declared as observed more than once (' varobs{i} ')!'])
        end
    end
end