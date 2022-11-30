@#include "backward_model.inc"

@#define simlength = 100

initialconditions = dseries(rand(2, @{5*n}), '2021Q3', { ...
@#for i in 1:n
  'x@{i}', 'y@{i}', 'z@{i}', 't@{i}', 's@{i}', ...
@#endfor
});

innovations = dseries(randn(@{simlength}, @{5*n}), '2022Q1', { ...
@#for i in 1:n
  'e_x@{i}', 'e_y@{i}', 'e_z@{i}', 'e_t@{i}', 'e_s@{i}', ...
@#endfor
});
