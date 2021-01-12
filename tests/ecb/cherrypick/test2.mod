var x y z;

model;

x = y(-1)+z(1);

y = x-z;

diff(z) = 0;

end;


if ~isfield(M_, 'exo_names')
   error('M_.exo_names not defined.')
end

if ~iscell(M_.exo_names)
   error('M_.exo_names not a cell array.')
end

if ~isempty(M_.exo_names)
   error('M_.exo_names not an empty cell array.')
end

if ~isfield(M_, 'param_names')
   error('M_.param_names not defined.')
end

if ~iscell(M_.param_names)
   error('M_.param_names not a cell array.')
end

if ~isempty(M_.param_names)
   error('M_.param_names not an empty cell array.')
end