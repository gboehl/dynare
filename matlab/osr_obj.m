function [loss,vx,info,exit_flag]=osr_obj(x,i_params,i_var,weights)
% objective function for optimal simple rules (OSR). Deprecated
% interface. New one: osr_obj_1.m
%
% INPUTS
%   x                         vector           values of the parameters
%                                              over which to optimize
%   i_params                  vector           index of optimizing parameters in M_.params
%   i_var                     vector           variables indices
%   weights                   vector           weights in the OSRs
%
% OUTPUTS
%   loss                      scalar           loss function returned to solver
%   vx                        vector           variances of the endogenous variables
%   info                      vector           info vector returned by resol
%   exit_flag                 scalar           exit flag returned to solver
%
% SPECIAL REQUIREMENTS
%   none
% Copyright (C) 2005-2013 Dynare Team
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

[loss,info,exit_flag,vx,junk]=osr_obj_1(x,i_params,i_var,weights);