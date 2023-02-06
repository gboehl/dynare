var x y;

varexo ex ey;

model;
  y = x*x/2 + ey;
  x = ex;
end;

initval;
  y = 0;
  x = 0;
  ex = 0;
  ey = 0;
end;

steady;

shocks;
    var ex;
    periods 1 2 3 4;
    values 1 .5 .25 .125;
end;

perfect_foresight_setup(periods=4);
perfect_foresight_solver;
