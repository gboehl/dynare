function alpha = a2alpha(a)

% Computes the m alpha coefficients from the m a coefficients of the PAC model.
%
% INPUTS 
% - a      [double]   m×1 vector of coefficients.
%
% OUTPUTS 
% - alpha  [double]   m×1 vector of coefficients.
%
% NOTES 
%
%  Given the current estimate of the PAC parameters a₀, a₁, ..., aₘ₋₁, the routine does the following:
%
%  αₘ = aₘ₋₁
%  αₘ₋₁ = aₘ₋₂ - aₘ₋₁
%  ⋮
%  αᵢ = aᵢ₋₁ - aᵢ ≡ - Δaᵢ
%  ⋮
%  α₂ = a₁ - a₂
%  α₁ = a₀ - a₁ - 1
%
% Computed coefficients αᵢ define lag polynomial A(L) = 1 + α₁L + α₂L² + … + αₘ Lᵐ such that a₀ = A(1).

% Copyright © 2018-2024 Dynare Team
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

% Return an error if the input is not a vector
if ~isvector(a)
    error('Input argument has to be a vector of doubles!')
end 

% Get the number of PAC parameters (without the discount factor)
m = length(a);

% Initialize the vector of transformed PAC parameters.
alpha = zeros(m, 1);

% Compute the transformed parameters
alpha(m) = a(m);
alpha(2:m-1) = a(2:m-1)-a(3:m);
alpha(1) = a(1)-a(2)-1;
