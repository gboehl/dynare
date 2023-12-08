function dyn_latex_table(M_, options_, title, LaTeXtitle, headers, labels, values, label_width, val_width, val_precis, optional_header)
%function dyn_latex_table(M_, options_, title, LaTeXtitle, headers, labels, values, label_width, val_width, val_precis, optional_header)

% Copyright © 2015-2020 Dynare Team
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

if options_.noprint
    return
end

if length(headers) < 2
    error('headers length must be >= 2')
end

OutputDirectoryName = CheckPath('latex',M_.dname);

%delete dollars in label as they will be added automatically below
begin_dollar = cellfun (@(x)startsWith(x,'$'),labels,'UniformOutput',1);
end_dollar = cellfun (@(x)endsWith(x,'$'),labels,'UniformOutput',1);

if all(begin_dollar) && all(end_dollar)
    labels = cellfun(@(x)delete_dollar(x,'begin'),labels,'UniformOutput',0);
    labels = cellfun(@(x)delete_dollar(x,'end'),labels,'UniformOutput',0);
end


%delete dollars in headers as they will be added automatically below
begin_dollar = cellfun (@(x)startsWith(x,'$'),headers,'UniformOutput',1);
end_dollar = cellfun (@(x)endsWith(x,'$'),headers,'UniformOutput',1);

if all(begin_dollar) && all(end_dollar)
    headers = cellfun(@(x)delete_dollar(x,'begin'),headers,'UniformOutput',0);
    headers = cellfun(@(x)delete_dollar(x,'end'),headers,'UniformOutput',0);
end

% Set width of label column
if isempty(label_width)
    label_width = cellofchararraymaxlength(vertcat(headers{1}, labels))+2;
else
    label_width = max(cellofchararraymaxlength(vertcat(headers{1}, labels))+2, label_width);
end
label_format_leftbound = sprintf('$%%-%ds$', label_width);

% Set width of other columns
if all(~isfinite(values))
    values_length = 4;
else
    values_length = max(ceil(max(max(log10(abs(values(isfinite(values))))))),1)+val_precis+1;
end
if any(values < 0) %add one character for minus sign
    values_length = values_length+1;
end
headers_length = cellofchararraymaxlength(headers(2:end));
if isempty(val_width)
    val_width = max(headers_length, values_length)+4;
else
    val_width = max(max(headers_length, values_length)+4, val_width);
end
value_format = sprintf('%%%d.%df', val_width, val_precis);
header_string_format = sprintf('$%%%ds$', val_width);

% Create and print header string
header_string = sprintf(label_format_leftbound, strrep(headers{1}, '\', '\\'));
header_code_string = ['l' repmat('c', 1, length(headers)-1)];
for i=2:length(headers)
    header_string = [header_string '\t & \t ' sprintf(header_string_format, strrep(headers{i},'\','\\'))];
end
header_string = [header_string '\\\\\n'];

filename = [OutputDirectoryName '/' M_.fname '_' LaTeXtitle '.tex'];
fidTeX = fopen(filename,'w');

stack = dbstack;
fprintf(fidTeX, ['%% ' datestr(now,0) ', created by ' stack(2).file]);
fprintf(fidTeX, ' \n');
fprintf(fidTeX, ' \n');
fprintf(fidTeX, '\\begin{center}\n');
fprintf(fidTeX, '\\begin{longtable}{%s} \n', header_code_string);
fprintf(fidTeX, ['\\caption{',title,'}\\\\\n ']);

fprintf(fidTeX, ['\\label{Table:',LaTeXtitle,'}\\\\\n']);
fprintf(fidTeX, '\\toprule \n');
if nargin==11
    for ii = 1:length(optional_header)
        fprintf(fidTeX,'%s\n',optional_header{ii});
    end
end
fprintf(fidTeX, header_string);
fprintf(fidTeX, '\\midrule \\endfirsthead \n');
fprintf(fidTeX, '\\caption{(continued)}\\\\\n ');
fprintf(fidTeX, '\\toprule \\\\ \n');
if nargin==11
    for ii = 1:length(optional_header)
        fprintf(fidTeX, '%s\n', optional_header{ii});
    end
end
fprintf(fidTeX, header_string);
fprintf(fidTeX, '\\midrule \\endhead \n');
fprintf(fidTeX, ['\\midrule \\multicolumn{',num2str(size(headers,1)),'}{r}{(Continued on next page)} \\\\ \\bottomrule \\endfoot \n']);
fprintf(fidTeX, '\\bottomrule \\endlastfoot \n');
for i = 1:size(values,1)
    fprintf(fidTeX, label_format_leftbound, labels{i});
    fprintf(fidTeX, ['\t & \t' value_format], values(i,:));
    fprintf(fidTeX, ' \\\\ \n');
end

fprintf(fidTeX, '\\end{longtable}\n ');
fprintf(fidTeX, '\\end{center}\n');
fprintf(fidTeX, '%% End of TeX file.\n');
fclose(fidTeX);

function x=delete_dollar(x,position_string)
if strcmp(position_string,'begin')
    x(1)=[];
elseif strcmp(position_string,'end')
    x(end)=[];
end

