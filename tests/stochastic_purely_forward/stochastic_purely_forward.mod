var a;
varexo epsil ;
parameters betta;
betta = 0.97;

model;
a = betta*a(+1)+epsil;
end;

initval;
a=0;
end;
steady;

shocks;
var epsil; stderr .2;
end;

steady;
check;

% Skip test under Octave 5.1
% ordeig() is buggy in that version (but is fixed in later ones; and in older
% ones it is absent, so we use our replacement)
if ~isoctave || octave_ver_less_than('5.1') || ~octave_ver_less_than('5.2')

stoch_simul(periods=0, irf=30, order=1);
stoch_simul(periods=2000, irf=30, order=1);

stoch_simul(periods=0, irf=30, order=1,hp_filter=1600);
stoch_simul(periods=2000, irf=30, order=1,hp_filter=1600);

end
