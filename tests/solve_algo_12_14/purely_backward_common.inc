@#include "backward_model.inc"

initval;
@#for i in 1:n
y@{i} = 1;
z@{i} = 0.5;
@#endfor
end;

shocks;
@#for i in 1:n
var e_x@{i};
periods 1;
values @{0.1/n};

var e_y@{i};
periods 2;
values @{0.05+0.1/n};

var e_z@{i};
periods 3;
values @{0.1+0.1/n};

var e_t@{i};
periods 4;
values @{0.15+0.1/n};

var e_s@{i};
periods 5;
values @{0.2+0.1/n};
@#endfor
end;

perfect_foresight_setup(periods = 100);
