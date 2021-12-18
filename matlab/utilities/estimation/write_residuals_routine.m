function write_residuals_routine(lhs, rhs, eqname, ipnames, DynareModel, pacmodl)

% Creates a routine for evaluating the residuals of the nonlinear equation.

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

fun = sprintf('r_%s', eqname);
fid = fopen(['+' DynareModel.fname filesep() fun '.m'], 'w');
fprintf(fid, 'function r = %s(params, data, DynareModel, DynareOutput)\n', fun);
fprintf(fid, '\n');
fprintf(fid, '%% Evaluates the residuals for equation %s.\n', eqname);
fprintf(fid, '%% File created by Dynare (%s).\n', datetime);
fprintf(fid, '\n');
for i=1:length(ipnames)
    fprintf(fid, 'DynareModel.params(%u) = params(%u);\n', ipnames(i), i);
end
fprintf(fid, '\n');
if nargin>5
    fprintf(fid, 'DynareModel = pac.update.parameters(''%s'', DynareModel, DynareOutput, false);\n', pacmodl);
    fprintf(fid, '\n');
end
fprintf(fid, 'r = %s-(%s);\n', lhs, rhs);
fclose(fid);