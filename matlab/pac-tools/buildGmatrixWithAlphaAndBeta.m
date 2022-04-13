function [G, alpha, beta] = buildGmatrixWithAlphaAndBeta(params)
    
% Builds the G matrix needed for PAC.
%
% INPUTS 
% - params    [double]    (m+1)*1 vector of PAC parameters.
%
% OUTPUTS 
% - G         [double]    (m+1)*(m+1) matrix.
% - alpha     [double]    m*1 vector of PAC parameters.
% - beta      [double]    scalar, discount factor.

% Copyright © 2018 Dynare Team
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

% Return an error if the input is not a vector.
if ~isvector(params) || ~isnumeric(params) || ~isreal(params)
    error('Input argument has to be a vector of doubles!')
end 

% Get the number of parameters
m = length(params)-1;

% Return an error if params is too small.
if m<1
    error('Input argument has to be a vector with at least two elements. The last element is the discount factor.')
end

% Get the transformed PAC parameters and discount factor.
alpha = flip(a2alpha(params(1:m)));
beta = params(end);

% Return an error if beta is not a discount factor
if beta<eps || beta>1-eps
    error('beta has to be a discount factor!')
end 

% Initialize the returned G matrix.
G = zeros(m);

% Fill the returned G matrix.
G(1:m-1,2:m) = eye(m-1);
G(m, :) = -transpose(alpha);
G(m, :) = G(m, :).*flip(cumprod(beta*ones(1,m)));