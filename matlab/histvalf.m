function [endo_histval, exo_histval, exo_det_histval] = histvalf(M, options)
%function [endo_histval, exo_histval, exo_det_histval] = histvalf(M, options)
% Sets initial values for simulation using values contained in `fname`, a
% file possibly created by a call to `smoother2histval`
%
% INPUTS
%    fname:                       name of file containing initial values
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none


% Copyright (C) 2014-2020 Dynare Team
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

series = histvalf_initvalf('HISTVAL', M, options);
% capture the difference between stochastic and
% perfect foresight setup
k = M.orig_maximum_lag - M.maximum_lag + 1;
endo_histval  = series{M.endo_names{:}}.data(k:end, :)';

exo_histval  = [];
if M.exo_nbr
    exo_histval  = series{M.exo_names{:}}.data(k:end, :)';
end
exo_det_histval  = [];
if M.exo_det_nbr
    exo_det_histval  = series{M.exo_names{:}}.data(k:end, :)';
end

