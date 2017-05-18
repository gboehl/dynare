function [nodes,weights] = spherical_radial_sigma_points(n)
%
% Computes nodes and weigths from a third-degree spherical-radial cubature
% rule.
% INPUTS
%    n                  [integer]   scalar, number of variables.
%
% OUTPUTS
%    nodes          [double]    nodes of the cubature
%    weigths        [double]    associated weigths
%
% REFERENCES
%
% Arasaratnam & Haykin 2008,2009.
%
% NOTES
%
% Copyright (C) 2009-2017 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

nodes = (sqrt(n)*([eye(n) -eye(n)]))' ;
weights = (1/(2*n)) ;
