function tf = ver_less_than(ver1, ver2)
%function tf = ver_less_than(ver1, ver2)
% ver1 < ver2 ? 1 : 0;
%
% INPUTS
%    ver1    [string]    software version number
%    ver2    [string]    software version number
%
% OUTPUTS
%    tf      [bool]      true if ver1 < ver2
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2015-2022 Dynare Team
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

if strcmp(ver1,ver2)
    tf = false;
    return
else
    tf = true;
end
ver1 = strsplit(ver1, {'.', '-'});
ver2 = strsplit(ver2, {'.', '-'});

maj_ver1 = str2double(ver1{1});
maj_ver2 = str2double(ver2{1});
if maj_ver1 < maj_ver2
    tf = true;
    return
elseif maj_ver1 > maj_ver2
    tf = false;
    return
end

min_ver1 = str2double(ver1{2});
min_ver2 = str2double(ver2{2});
if (maj_ver1 == maj_ver2) && (min_ver1 < min_ver2)
    tf = true;
    return
elseif (maj_ver1 == maj_ver2) && (min_ver1 > min_ver2)
    tf = false;
    return
end

%deal with revision in Dynare 4 and unstable versions involved
if min(maj_ver1,maj_ver2)<5 %old versioning scheme with three digits
    if (length(ver1) == length(ver2) && length(ver1) == 3)
        %check if master branch (unstable) or stable
        ismaster1 = isnan(str2double(ver1{3}));
        ismaster2 = isnan(str2double(ver2{3}));
        if (maj_ver1 == maj_ver2) && (min_ver1 == min_ver2) && (~ismaster1 && ismaster2)
            %ver2 is the unstable
            return
        end

        if ~ismaster1 && ~ismaster2 %both are stable versions
            rev_ver1 = str2double(ver1{3});
            rev_ver2 = str2double(ver2{3});
            if (maj_ver1 == maj_ver2) && (min_ver1 == min_ver2) && (rev_ver1 < rev_ver2)
                %ver1 has the lower minor version
                return
            end
        end
    else
        %ver1 is an unstable version
        error('Case is undefined, please contact the developers')
    end
elseif min(maj_ver1,maj_ver2)>=5 %new versioning scheme with three digits
    if strcmp(ver1{2},'x') || strcmp(ver1{2},'unstable')
        date_number_1=datenum([ver1{3} ver1{4} ver1{5}],'YYYYMMDD');
        stable_version_indicator_1=0;
    elseif ~isnan(str2double(ver1{2}))
        stable_version_indicator_1=1;
    else
        error('Case is undefined, please contact the developers')
    end
    if strcmp(ver2{2},'x') || strcmp(ver2{2},'unstable')
        date_number_2=datenum([ver2{3} ver2{4} ver2{5}],'YYYYMMDD');
        stable_version_indicator_2=0;
    elseif ~isnan(str2double(ver2{2}))
        stable_version_indicator_2=1;
    end
    if ~stable_version_indicator_1 && ~stable_version_indicator_2
        if date_number_1<date_number_2
            return
        end
    else
        %comparison between unstable and stable version
        error('You cannot compare a stable release to an unstable version of the same branch.')
    end
end
tf = false;

return % --*-- Unit tests --*--

%@test:1
ver1='4.4';
ver2='4.5.2';
t(1)=dassert(ver_less_than(ver1,ver2),true);
T = all(t);
%@eof:1

%@test:2
ver1='4.4';
ver2='6-unstable-2021-12-15-1737-21a8a579';
t(1)=dassert(ver_less_than(ver1,ver2),true);
T = all(t);
%@eof:2

%@test:3
ver1='5.0';
ver2='5.1';
t(1)=dassert(ver_less_than(ver1,ver2),true);
T = all(t);
%@eof:3

%@test:4
ver1='6-unstable-2021-12-18-1227-c43777f6';
ver2='6-unstable-2021-12-19-1953-d841fc7c';
t(1)=dassert(ver_less_than(ver1,ver2),true);
T = all(t);
%@eof:4

%@test:5
ver1='5.5';
ver2='5.5';
t(1)=dassert(ver_less_than(ver1,ver2),false);
T = all(t);
%@eof:5
