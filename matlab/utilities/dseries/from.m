function from(varargin)

% Copyright (C) 2014 Dynare Team
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

lvarargin = lower(varargin);

if ~(ismember('to',lvarargin) && ismember('do',lvarargin))
    error(get_error_message_0())
end

to_id = strmatch('to',lvarargin);
do_id = strmatch('do',lvarargin);

if do_id<to_id
    msg = sprinf('Wrong syntax! The TO keyword must preceed the DO keyword.\n');
    error(get_error_message_0(msg))
end

if ~isdate(varargin{1})
    msg = sprintf('Wrong syntax! The FROM statement must be followed by a dates object.\n');
    error(get_error_message_0(msg))
end

if ~isequal(to_id,2)
    msg = sprintf('Wrong syntax! The first dates object must be immediately followed by the TO keyword.\n');
    error(get_error_message_0(msg))
end

if ~isdate(varargin{3})
    msg = sprintf('Wrong syntax! The TO keyword must be followed by a second dates object.\n');
    error(get_error_message_0(msg))
end

d1 = dates(varargin{1});
d2 = dates(varargin{3});

if d1>d2
    error('The first date must preceed the second one!')
end

if ~isequal(do_id,4)
    msg = sprintf('Wrong syntax! The second dates object must be immediately followed by the DO keyword.\n');
    error(get_error_message_0(msg))
end

% Build the recursive expression.
EXPRESSION = [];
for i=5:nargin
    EXPRESSION = [EXPRESSION, varargin{i}];
end

% Get all the variables involved in the recursive expression.
variables = regexpi(EXPRESSION, '\w*\(t\)|\w*\(t\-\d\)|\w*\(t\+\d\)','match');

% Remove the time indices.
for i=1:length(variables)
    start = regexpi(variables{i},'\(t\)|\(t\-\d\)|\(t\+\d\)');
    variables(i) = {variables{i}(1:start-1)};
end

% Remove duplicates.
variables = unique(variables);

% Test that all the involved variables are available dseries objects.
for i=1:length(variables)
    try
        var = evalin('caller',variables{i});
    catch
        error(['Variable ' variables{i} ' is unknown!'])
    end
    if ~isdseries(var)
        error(['Variable ' variables{i} ' is not a dseries object!'])
    else
        eval(sprintf('%s = var;',variables{i}));
    end
end

% Check that the recursion is assigning something to a variable
equal_id = strfind(EXPRESSION,'=');
if isempty(equal_id)
    error('The expression following the DO keyword must be an assignment (missing equal symbol)!')
end
if isequal(length(equal_id),1)
    assignedvariablename = regexpi(EXPRESSION(1:equal_id-1), '\w*\(t\)|\w*\(t\-\d\)|\w*\(t\+\d\)','match');
    if isempty(assignedvariablename)
        error('The expression following the DO keyword must be an assignment (missing variable before the equal symbol)!')
    end
    if length(assignedvariablename)>1
        error('No more than one variable can be assigned!')
    end
    assignedvariablename = assignedvariablename{1}(1:regexpi(assignedvariablename{1},'\(t\)|\(t\-\d\)|\(t\+\d\)')-1);
    eval(sprintf('wrongtype = ~isdseries(%s);',assignedvariablename))
    if wrongtype
        error('The assigned variable must be a dseries object!')
    end
else
    error('Not yet implemented! Only one assignment is allowed in the FROM-TO-DO statement.')
end

% Transform the indexed variables after the assignment symbol: X(t-1) -> X(t-1).data
expression = EXPRESSION(equal_id+1:end);
expression = regexprep(expression,'\w*\(t\)|\w*\(t\-\d\)|\w*\(t\+\d\)','$0.data');
EXPRESSION = [EXPRESSION(1:equal_id),expression];

% Run the recursion!
eval(sprintf('t=dates(''%s''); while t<=dates(''%s''), %s; t = t+1; end',char(d1),char(d2),EXPRESSION))

% Put assigned variable back in the caller workspace...
eval(sprintf('assignin(''caller'', assignedvariablename, %s)',assignedvariablename));



function msg = get_error_message_0(msg)
    if ~nargin
        msg = sprintf('Wrong syntax! The correct syntax is:\n\n');
    else
        msg = [msg, sprintf('The correct syntax is:\n\n')];
    end
    msg = [msg, sprintf('    from d1 to d2 do SOMETHING\n\n')];
    msg = [msg, sprintf('where d1<d2 are dates objects, and SOMETHING is a recursive expression involving dseries objects.')];
