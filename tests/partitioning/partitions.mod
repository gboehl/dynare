var y (country='US', long_name='United States') z (country='EA', long_name='Euro Area') x (long_name='March');

varexo e (country='US', kind='toto')
       f (country='US', kind='titi')
       g (country='EA', kind='toto')
       h (country='EA', kind='titi') ;

parameters aUS (country='US', colour='red')
           aEA (country='EA', colour='red')
           bUS (country='US', colour='blue')
           bEA (country='EA', colour='blue')
           cUS (country='US', colour='yellow')
           cEA (country='EA', colour='yellow') ;

model;
y = .5*y(-1) + e-f^2 ;
z = .5*z(-1) + g-h^2 ;
x = e-g;
end;

if ~(isfield(M_.exo_partitions, 'country') && isequal(M_.exo_partitions.country, {'US'  'US'  'EA'  'EA'}))
    error('First partition on exogenous variables is wrong!')
end

if ~(isfield(M_.exo_partitions, 'kind') && isequal(M_.exo_partitions.kind, {'toto'  'titi'  'toto'  'titi'}))
    error('Second partition on exogenous variables is wrong!')
end

if ~(isfield(M_.endo_partitions, 'country') && isequal(M_.endo_partitions.country, {'US'  'EA' ''}))
    error('Partition on endogenous variables is wrong!')
end

if ~(isfield(M_.param_partitions, 'colour') && isequal(M_.param_partitions.colour, {'red'  'red' 'blue'  'blue' 'yellow'  'yellow'}))
    error('First partition on parameters is wrong!')
end

if ~(isfield(M_.param_partitions, 'country') && isequal(M_.param_partitions.country, {'US'  'EA' 'US'  'EA' 'US'  'EA'}))
    error('Second partition on parameters is wrong!')
end

if isfield(M_.endo_partitions, 'long_name')
    error('long_name should not be a member of endo_partitions!')
end

tmp = {'United States'; 'Euro Area'; 'March'};

if ~isequal(M_.endo_names_long, tmp)
    error('M_.endo_names_long is wrong!')
end
