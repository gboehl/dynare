function write_latex_definitions
%function write_latex_definitions
% Writes a latex file containing the variable names, latex names, and
% tags/comments
%
% INPUTS
%    none
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2013 Dynare Team
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

global M_;

if M_.exo_det_nbr == 0
    tables = {'Endogenous', 'Exogenous', 'Parameters'};
    M_var_root = {'M_.endo', 'M_.exo', 'M_.param'};
else
    tables = {'Endogenous', 'Exogenous', 'Exogenous Deterministic', 'Parameters'};
    M_var_root = {'M_.endo', 'M_.exo', 'M_.exo_det', 'M_.param'};
end
fid = fopen([M_.fname '_latex_definitions.tex'], 'w');
fprintf(fid, '\\documentclass[10pt,a4paper]{article}\n');
fprintf(fid, '\\usepackage{geometry}\n');
fprintf(fid, '\\begin{document}\n');

for i=1:length(tables)
    fprintf(fid, '\\begin{table}[ht]\n');
    fprintf(fid, ['\\caption{' tables{i} '}\n']);
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\begin{tabular}{c c c}\n');
    fprintf(fid, '\\hline\\hline\n');
    fprintf(fid, 'Variable & LaTeX & Description\\\\\n');
    fprintf(fid, '\\hline\n');

    names = eval([M_var_root{i} '_names']);
    tex = eval([M_var_root{i} '_names_tex']);
    long = eval([M_var_root{i} '_names_long']);
    for j=1:size(names,1)
        fprintf(fid, '%s & $%s$ & %s\\\\\n', ...
            regexprep(strtrim(names(j,:)), '_', '\\_'), ...
            strtrim(tex(j,:)), ...
            regexprep(strtrim(long(j,:)), '_', '\\_'));
    end

    fprintf(fid, '\\hline\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');
    fprintf(fid, '\\newpage\n');
end

fprintf(fid, '\\end{document}\n');
fclose(fid);
end
