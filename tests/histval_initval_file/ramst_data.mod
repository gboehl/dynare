/* Verify that the “datafile” option of “perfect_foresight_setup” behaves as
   “initval_file” (see #1663) */

var c k;
varexo x;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;


model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
end;

initval;
x = 1;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
end;

steady;

shocks;
  var x;
  periods 2;
  values 1.1;
end;

perfect_foresight_setup(periods=200);
perfect_foresight_solver;

fh = fopen('ramst_data.m', 'w');
fprintf(fh, 'INIT__ = ''1Y'';\n');
fprintf(fh, 'NAMES__ = {''c'', ''k'', ''x''};\n');
fprintf(fh, 'TEX__ = {''c'', ''k'', ''x''};\n');
fprintf(fh, 'c = [');
fprintf(fh, '%f ', oo_.endo_simul(1,:));
fprintf(fh, '];\n');
fprintf(fh, 'k = [');
fprintf(fh, '%f ', oo_.endo_simul(2,:));
fprintf(fh, '];\n');
fprintf(fh, 'x = [');
fprintf(fh, '%f ', oo_.exo_simul);
fprintf(fh, '];\n');
fclose(fh);

INIT__ = '1Y';
NAMES__ = {'c', 'k', 'x'};
TEX__  = {'c', 'k', 'x'};
eval('c = oo_.endo_simul(1,:);');
eval('k = oo_.endo_simul(2,:);');
eval('x = oo_.exo_simul'';');
save('ramst_data.mat', 'INIT__', 'NAMES__', ...
     'TEX__', 'c', 'k', 'x');

fh = fopen('ramst_data.csv', 'w');
fprintf(fh, 'c,k,x\n');
for i = 1:size(oo_.endo_simul, 2);
  fprintf(fh, '%f, ', oo_.endo_simul(:, i));
  fprintf(fh, '%f\n', oo_.exo_simul(i));
end;
fclose(fh);

