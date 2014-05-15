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

to_id = strmatch('to',lvarargin);
do_id = strmatch('do',lvarargin);

if isempty(to_id) || isempty(do_id)
    error(get_error_message_0())
end

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

d1 = dates(varargin{1}); % First date
d2 = dates(varargin{3}); % Last date

if d1>d2
    error('The first date must preceed the second one!')
end

if ~isequal(do_id,4)
    msg = sprintf('Wrong syntax! The second dates object must be immediately followed by the DO keyword.\n');
    error(get_error_message_0(msg))
end

% Build the recursive expression.
EXPRESSION = char([varargin{5:end}]);

% Get all the variables involved in the recursive expression.
variables = unique(regexpi(EXPRESSION, '\w*\(t\)|\w*\(t\-\d\)|\w*\(t\+\d\)','match'));

% Build an incidence table (max lag/lead for each variable)
%
% Column 1: Name of the variable.
% Column 2: Maximum lag order.
% Column 3: Equal to 1 if the variable appears at the current period, 0 otherwise.
% Column 4: Maximum lead order.
% Column 5: Vector of effective lag orders.
% Column 6: Vector of effective lead orders.
%
% Initialization.
leadlagtable = cell(0,6);
% Loop over the variables (dseries objects).
for i=1:length(variables)
    current = ~isempty(regexpi(variables{i},'\(t\)'));
    lag = ~isempty(regexpi(variables{i},'\(t\-\d\)'));
    lead = ~isempty(regexpi(variables{i},'\(t\+\d\)'));
    start = regexpi(variables{i},'\(t\)|\(t\-\d\)|\(t\+\d\)');
    index = variables{i}(start:end);
    variables(i) = {variables{i}(1:start-1)};
    if isempty(leadlagtable)
        leadlagtable(1,1) = {variables{i}};
        if current
            leadlagtable(1,3) = {1};
        else
            leadlagtable(1,3) = {0};
        end
        if lag
            tmp = regexpi(index,'\d','match');
            leadlagtable(1,2) = {str2double(tmp{1})};
            leadlagtable(1,5) = {str2double(tmp{1})};
        else
            leadlagtable(1,2) = {0};
            leadlagtable(1,5) = {[]};
        end
        if lead
            tmp = regexpi(index,'\d','match');
            leadlagtable(1,4) = {str2double(tmp{1})};
            leadlagtable(1,6) = {str2double(tmp{1})};
        else
            leadlagtable(1,4) = {0};
            leadlagtable(1,6) = {[]};
        end
    else
        linea = strmatch(variables{i},leadlagtable(:,1));
        if isempty(linea)
            % This is a new variable!
            linea = size(leadlagtable,1)+1;
            leadlagtable(linea,1) = {variables{i}};
            leadlagtable(linea,2) = {0};
            leadlagtable(linea,3) = {0};
            leadlagtable(linea,4) = {0};
            leadlagtable(linea,5) = {[]};
            leadlagtable(linea,6) = {[]};
        end
        if current
            leadlagtable(linea,3) = {1};
        end
        if lag
            tmp = regexpi(index,'\d','match');
            leadlagtable(linea,2) = {max(str2double(tmp{1}),leadlagtable{linea,2})};
            leadlagtable(linea,5) = {sortrows([leadlagtable{linea,5}; str2double(tmp{1})])};
        end
        if lead
            tmp = regexpi(index,'\d','match');
            leadlagtable(linea,4) = {max(str2double(tmp{1}),leadlagtable{linea,4})};
            leadlagtable(linea,6) = {sortrows([leadlagtable{linea,6}; str2double(tmp{1})])};
        end
    end
end

% Set the number of variables
number_of_variables = size(leadlagtable,1);

% Test that all the involved variables are available dseries objects. Also check that
% these time series are defined over the time range given by d1 and d2 (taking care of
% the lags and leads).
for i=1:number_of_variables
    current_variable = leadlagtable{i,1};
    try
        var = evalin('caller',current_variable);
    catch
        error(['dseries::from: Variable ' current_variable ' is unknown!'])
    end
    if ~isdseries(var)
        error(['dseries::from: Variable ' current_variable ' is not a dseries object!'])
    else
        if d1<var.dates(1)+leadlagtable{i,2}
            msg = sprintf('dseries::from: Initial date of the loop (%s) is inconsistent with %s''s range!\n',char(d1),current_variable);
            msg = [msg, sprintf('               Initial date should be greater than or equal to %s.',char(var.dates(1)+leadlagtable{i,2}))];
            error(msg)
        end
        if d2>var.dates(end)-leadlagtable{i,4}
            msg = sprintf('dseries::from: Last date of the loop (%s) is inconsistent with %s''s range!\n',char(d2),current_variable);
            msg = [msg, sprintf('               Last date should be less than or equal to %s.',char(var.dates(end)-leadlagtable{i,4}))];
            error(msg)
        end
        eval(sprintf('%s = var;',current_variable));
    end
end

% Check that the recursion is assigning something to a variable
equal_id = strfind(EXPRESSION,'=');
if isempty(equal_id)
    error('Wrong syntax! The expression following the DO keyword must be an assignment (missing equal symbol).')
end
if isequal(length(equal_id),1)
    % Get the name of the assigned variable (with time index)
    assignedvariablename = regexpi(EXPRESSION(1:equal_id-1), '\w*\(t\)|\w*\(t\-\d\)|\w*\(t\+\d\)','match');
    if isempty(assignedvariablename)
        error('Wrong syntax! The expression following the DO keyword must be an assignment (missing variable before the equal symbol).')
    end
    if length(assignedvariablename)>1
        error('No more than one variable can be assigned!')
    end
    % Check that the dynamic model for the endogenous variable is not forward looking.
    start = regexpi(assignedvariablename{1},'\(t\)|\(t\-\d\)|\(t\+\d\)');
    index = assignedvariablename{1}(start:end);
    assignedvariablename = assignedvariablename{1}(1:start-1);
    indum = index2num(index);
    indva = strmatch(assignedvariablename, leadlagtable(:,1));
    if indum<leadlagtable{indva,4}
        error('dseries::from: It is not possible to simulate a forward looking model!')
    end
    % Check that the assigned variable does not depend on itself (the assigned variable can depend on its past level but not on the current level).
    tmp = regexpi(EXPRESSION(equal_id+1:end), sprintf('%s\\(t\\)|%s\\(t\\-\\d\\)|%s\\(t\\+\\d\\)',assignedvariablename,assignedvariablename,assignedvariablename),'match');
    tmp = cellfun(@extractindex, tmp);
    tmp = cellfun(@index2num, tmp);
    if ~all(tmp(:)<indum)
        error(sprintf('dseries::from: On the righthand side, the endogenous variable, %s, must be indexed by %s at most.',assignedvariablename,num2index(indum-1)))
    end
else
    error('Not yet implemented! Only one assignment is allowed in the FROM-TO-DO statement.')
end

% Put all the variables in a unique dseries object.
list_of_variables = leadlagtable{1,1};
for i=2:number_of_variables
   list_of_variables = [list_of_variables, ',' leadlagtable{i,1}];
end
try
    eval(sprintf('tmp = [%s];', list_of_variables));
catch
    error('dseries::from: All the dseries objects should contain variables with different names!')
end

% Get base time index
t1 = find(d1==tmp.dates);
t2 = find(d2==tmp.dates);

% Get data
data = tmp.data;

% Transform EXPRESSION by replacing calls to the dseries objects by references to data.
for i=1:number_of_variables
    EXPRESSION = regexprep(EXPRESSION,sprintf('%s\\(t\\)',leadlagtable{i,1}),sprintf('data(t,%s)',num2str(i)));
    for j=1:length(leadlagtable{i,5})
        lag = leadlagtable{i,5}(j);
        EXPRESSION = regexprep(EXPRESSION,sprintf('%s\\(t-%s\\)',leadlagtable{i,1},num2str(lag)),sprintf('data(t-%s,%s)',num2str(lag),num2str(i)));
    end
    for j=1:length(leadlagtable{i,6})
        lead = leadlagtable{i,6}(j);
        EXPRESSION = regexprep(EXPRESSION,sprintf('%s\\(t+%s\\)',leadlagtable{i,1},num2str(lead)),sprintf('data(t+%s,%s)',num2str(lead),num2str(i)));
    end
end

% Do the job. Evaluate the recursion.
eval(sprintf('for t=%s:%s, %s; end',num2str(t1),num2str(t2),EXPRESSION));

% Put assigned variable back in the caller workspace...
eval(sprintf('assignin(''caller'', ''%s'', dseries(data(:,indva),y.init,y.name,y.tex));',assignedvariablename))



function msg = get_error_message_0(msg)
    if ~nargin
        msg = sprintf('Wrong syntax! The correct syntax is:\n\n');
    else
        msg = [msg, sprintf('The correct syntax is:\n\n')];
    end
    msg = [msg, sprintf('    from d1 to d2 do SOMETHING\n\n')];
    msg = [msg, sprintf('where d1<d2 are dates objects, and SOMETHING is a recursive expression involving dseries objects.')];


function index = extractindex(str)
    index = regexpi(str,'\(t\)|\(t\-\d\)|\(t\+\d\)','match');


function i = index2num(id)
    if isequal('(t)',id)
        i = 0;
        return
    end
    if isequal('-',id(3))
        i = - str2num(id(4:end-1));
    else
        i = str2num(id(4:end-1));
    end


function id = num2index(i)
    if isequal(i,0)
        id = '(t)';
        return
    end
    if i<0
        id = ['(t-' int2str(abs(i)) ')'];
    else
        id = ['(t+' int2str(i) ')'];
    end