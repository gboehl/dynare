function o = addData(o, varargin)
%function o = addData(o, varargin)
% Add data to the current section of the current page in the report
%
% INPUTS
%   o          [report]  report object
%   varargin             arguments to @report_table/addData
%
% OUTPUTS
%   o          [report]  updated report object
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2019 Dynare Team
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

assert(~isempty(o.pages) , ...
       ['@report.addData: Before adding data, you must add a page, ' ...
        'section, and a table.']);
assert(~isempty(o.pages{end}.sections) , ...
       ['@report.addData: Before adding data, you must add a section and ' ...
        'a table']);
assert(~isempty(o.pages{end}.sections{end}.elements), ...
       '@report.addData: Before adding data, you must add a table');
assert(isa(o.pages{end}.sections{end}.elements{end}, 'report_table'), ...
       '@report.addData: you can only add data to a report_table object');

o.pages{end}.sections{end}.elements{end} = ...
    o.pages{end}.sections{end}.elements{end}.addData(varargin{:});
end
