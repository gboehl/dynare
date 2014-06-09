function tex = name2tex(name, info) % --*-- Unitary tests --*--

% Converts plain text name into tex name.

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

if nargin<2
    info = 0;
end

if info
    if iscell(name)
        nn = length(name);
    else
        nn = 1;
    end
end

if isoctave
    tex = regexprep(name, '_', '\_');
else
    tex = regexprep(name, '_', '\\_');
end

if info
    for i=1:nn
        if iscell(name)
            texname = tex{i};
        else
            texname = tex;
        end
        idx = strfind(texname,'_');
        ndx = length(idx);
        ntx = length(texname);
        if ndx
            gotonextcondition = 1;
            if isequal(ndx,1) && ~isequal(idx,2) && ~isequal(idx,ntx)
                texname = [ texname(1:idx-2) '_{' texname(idx+1:end) '}'];
                gotonextcondition = 0;
            end
            if gotonextcondition && isequal(ndx,2) && ~isequal(idx(1),2) && isequal(idx(2),ntx)
                texname = [ texname(1:idx(1)-2) '_{' texname(idx(1)+1:end) '}' ];
                gotonextcondition = 0;
            end
            if gotonextcondition && isequal(ndx,2) && idx(2)<ntx
                texname = [ texname(1:idx(2)-2) '_{' texname(idx(2)+1:end) '}' ];
                gotonextcondition = 0;
            end
            if gotonextcondition && ndx>2
                if idx(end)<ntx
                    texname = [ texname(1:idx(end)-2) '_{' texname(idx(end)+1:end) '}' ];
                else
                    texname = [ texname(1:idx(end-1)-2) '_{' texname(idx(end-1)+1:end) '}' ];
                end
            end
            if iscell(name)
                tex(i) = { texname };
            else
                tex = texname;
            end
        end
    end
end

%@test:1
%$ t = zeros(16,1);
%$ t1 = name2tex('_azert');
%$ t2 = name2tex('azert_');
%$ t3 = name2tex('_azert_');
%$ t4 = name2tex('azert_uiop');
%$ t5 = name2tex('azert_uiop_qsdfg');
%$ t6 = name2tex('azert_uiop_qsdfg_');
%$ t7 = name2tex('_azert_uiop_qsdfg');
%$ t8 = name2tex('_azert_uiop_qsdfg_');
%$ t11 = name2tex('_azert',1);
%$ t12 = name2tex('azert_',1);
%$ t13 = name2tex('_azert_',1);
%$ t14 = name2tex('azert_uiop',1);
%$ t15 = name2tex('azert_uiop_qsdfg',1);
%$ t16 = name2tex('azert_uiop_qsdfg_',1);
%$ t17 = name2tex('_azert_uiop_qsdfg',1);
%$ t18 = name2tex('_azert_uiop_qsdfg_',1);
%$
%$ t(1) = dyn_assert(strcmp(t1,'\\_azert'),1);
%$ t(2) = dyn_assert(strcmp(t2,'azert\\_'),1);
%$ t(3) = dyn_assert(strcmp(t3,'\\_azert\\_'),1);
%$ t(4) = dyn_assert(strcmp(t4,'azert\\_uiop'),1);
%$ t(5) = dyn_assert(strcmp(t5,'azert\\_uiop\\_qsdfg'),1);
%$ t(6) = dyn_assert(strcmp(t6,'azert\\_uiop\\_qsdfg\\_'),1);
%$ t(7) = dyn_assert(strcmp(t7,'\\_azert\\_uiop\\_qsdfg'),1);
%$ t(8) = dyn_assert(strcmp(t8,'\\_azert\\_uiop\\_qsdfg\\_'),1);
%$ t(9) = dyn_assert(strcmp(t11,'\\_azert'),1);
%$ t(10) = dyn_assert(strcmp(t12,'azert\\_'),1);
%$ t(11) = dyn_assert(strcmp(t13,'\\_azert\\_'),1);
%$ t(12) = dyn_assert(strcmp(t14,'azert_{uiop}'),1);
%$ t(13) = dyn_assert(strcmp(t15,'azert\\_uiop_{qsdfg}'),1);
%$ t(14) = dyn_assert(strcmp(t16,'azert\\_uiop_{qsdfg\\_}'),1);
%$ t(15) = dyn_assert(strcmp(t17,'\\_azert\\_uiop_{qsdfg}'),1);
%$ t(16) = dyn_assert(strcmp(t18,'\\_azert\\_uiop_{qsdfg\\_}'),1);
%$
%$ T = all(t);
%@eof:1

%@test:2
%$ t = zeros(16,1);
%$ t1 = name2tex({'_azert'});
%$ t2 = name2tex({'azert_'});
%$ t3 = name2tex({'_azert_'});
%$ t4 = name2tex({'azert_uiop'});
%$ t5 = name2tex({'azert_uiop_qsdfg'});
%$ t6 = name2tex({'azert_uiop_qsdfg_'});
%$ t7 = name2tex({'_azert_uiop_qsdfg'});
%$ t8 = name2tex({'_azert_uiop_qsdfg_'});
%$ t11 = name2tex({'_azert'},1);
%$ t12 = name2tex({'azert_'},1);
%$ t13 = name2tex({'_azert_'},1);
%$ t14 = name2tex({'azert_uiop'},1);
%$ t15 = name2tex({'azert_uiop_qsdfg'},1);
%$ t16 = name2tex({'azert_uiop_qsdfg_'},1);
%$ t17 = name2tex({'_azert_uiop_qsdfg'},1);
%$ t18 = name2tex({'_azert_uiop_qsdfg_'},1);
%$
%$ t(1) = dyn_assert(t1,{'\\_azert'});
%$ t(2) = dyn_assert(t2,{'azert\\_'});
%$ t(3) = dyn_assert(t3,{'\\_azert\\_'});
%$ t(4) = dyn_assert(t4,{'azert\\_uiop'});
%$ t(5) = dyn_assert(t5,{'azert\\_uiop\\_qsdfg'});
%$ t(6) = dyn_assert(t6,{'azert\\_uiop\\_qsdfg\\_'});
%$ t(7) = dyn_assert(t7,{'\\_azert\\_uiop\\_qsdfg'});
%$ t(8) = dyn_assert(t8,{'\\_azert\\_uiop\\_qsdfg\\_'});
%$ t(9) = dyn_assert(t11,{'\\_azert'});
%$ t(10) = dyn_assert(t12,{'azert\\_'});
%$ t(11) = dyn_assert(t13,{'\\_azert\\_'});
%$ t(12) = dyn_assert(t14,{'azert_{uiop}'});
%$ t(13) = dyn_assert(t15,{'azert\\_uiop_{qsdfg}'});
%$ t(14) = dyn_assert(t16,{'azert\\_uiop_{qsdfg\\_}'});
%$ t(15) = dyn_assert(t17,{'\\_azert\\_uiop_{qsdfg}'});
%$ t(16) = dyn_assert(t18,{'\\_azert\\_uiop_{qsdfg\\_}'});
%$
%$ T = all(t);
%@eof:2

%@test:3
%$ t = zeros(4,1);
%$ try
%$     t1 = name2tex({'_azert';'azert_';'_azert_';'azert_uiop';'azert_uiop_qsdfg';'azert_uiop_qsdfg_'});
%$     t(1) = 1;
%$ catch
%$     % Nothing to do here.
%$ end
%$
%$ if t(1)
%$     try
%$         t2 = name2tex({'_azert';'azert_';'_azert_';'azert_uiop';'azert_uiop_qsdfg';'azert_uiop_qsdfg_'},1);
%$         t(2) = 1;
%$     catch
%$         % Nothing to do here.
%$     end
%$ end
%$
%$ if t(1)
%$     t(3) = dyn_assert(t1,{'\\_azert';'azert\\_';'\\_azert\\_';'azert\\_uiop';'azert\\_uiop\\_qsdfg';'azert\\_uiop\\_qsdfg\\_'});
%$ end
%$
%$ if t(2)
%$     t(4) = dyn_assert(t2,{'\\_azert';'azert\\_';'\\_azert\\_';'azert_{uiop}';'azert\\_uiop_{qsdfg}';'azert\\_uiop_{qsdfg\\_}'});
%$ end
%$
%$ T = all(t);
%@eof:3

%@test:4
%$ try
%$     db = dseries('csv/dd.csv');
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ T = all(t);
%@eof:4
