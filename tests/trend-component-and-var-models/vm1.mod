var y1 y2 y3 y4;

varexo e1 e2 e3 e4;

parameters A111 A112 A113 A114
           A121 A122 A123 A124
           A131 A132 A133 A134
           A141 A142 A143 A144
           A211 A212 A213 A214
           A221 A222 A223 A224
           A231 A232 A233 A234
           A241 A242 A243 A244
           A311 A312 A313 A314
           A321 A322 A323 A324
           A331 A332 A333 A334
           A341 A342 A343 A344 ;

A1 = randn(4);
A2 = randn(4);
A3 = randn(4);

A111 = A1(1,1);
A112 = A1(1,2);
A113 = A1(1,3);
A114 = A1(1,4);
A121 = A1(2,1);
A122 = A1(2,2);
A123 = A1(2,3);
A124 = A1(2,4);
A131 = A1(3,1);
A132 = A1(3,2);
A133 = A1(3,3);
A134 = A1(3,4);
A141 = A1(4,1);
A142 = A1(4,2);
A143 = A1(4,3);
A144 = A1(4,4);

A211 = A2(1,1);
A212 = A2(1,2);
A213 = A2(1,3);
A214 = A2(1,4);
A221 = A2(2,1);
A222 = A2(2,2);
A223 = A2(2,3);
A224 = A2(2,4);
A231 = A2(3,1);
A232 = A2(3,2);
A233 = A2(3,3);
A234 = A2(3,4);
A241 = A2(4,1);
A242 = A2(4,2);
A243 = A2(4,3);
A244 = A2(4,4);

A311 = A3(1,1);
A312 = A3(1,2);
A313 = A3(1,3);
A314 = A3(1,4);
A321 = A3(2,1);
A322 = A3(2,2);
A323 = A3(2,3);
A324 = A3(2,4);
A331 = A3(3,1);
A332 = A3(3,2);
A333 = A3(3,3);
A334 = A3(3,4);
A341 = A3(4,1);
A342 = A3(4,2);
A343 = A3(4,3);
A344 = A3(4,4);

var_model(model_name=toto, eqtags=['eq:y1', 'eq:y2', 'eq:y3', 'eq:y4']);

model;

[name='eq:y1']
y1 = A111*y1(-1) + A112*y2(-1) + A113*y3(-1) + A114*y4(-1) +
     A211*y1(-2) + A212*y2(-2) + A213*y3(-2) + A214*y4(-2) +
     A311*y1(-3) + A312*y2(-3) + A313*y3(-3) + A314*y4(-3) + e1; 

[name='eq:y2']
y2 = A121*y1(-1) + A122*y2(-1) + A123*y3(-1) + A124*y4(-1) +
     A221*y1(-2) + A222*y2(-2) + A223*y3(-2) + A224*y4(-2) +
     A321*y1(-3) + A322*y2(-3) + A323*y3(-3) + A324*y4(-3) + e2; 

[name='eq:y3']
y3 = A131*y1(-1) + A132*y2(-1) + A133*y3(-1) + A134*y4(-1) +
     A231*y1(-2) + A232*y2(-2) + A233*y3(-2) + A234*y4(-2) +
     A331*y1(-3) + A332*y2(-3) + A333*y3(-3) + A334*y4(-3) + e3; 

[name='eq:y4']
y4 = A141*y1(-1) + A142*y2(-1) + A143*y3(-1) + A144*y4(-1) +
     A241*y1(-2) + A242*y2(-2) + A243*y3(-2) + A244*y4(-2) +
     A341*y1(-3) + A342*y2(-3) + A343*y3(-3) + A344*y4(-3) + e4;

end;

[A0, AR, B] = get_companion_matrix('toto');

if max(max(abs(AR(:,:,1)-A1)))>1e-12
   error('First order autoregressive matrix is wrong.')
end

if max(max(abs(AR(:,:,2)-A2)))>1e-12
   error('Second order autoregressive matrix is wrong.')
end

if max(max(abs(AR(:,:,3)-A3)))>1e-12
   error('Third order autoregressive matrix is wrong.')
end

CompanionMatrix = [A1, A2, A3;
                   eye(4), zeros(4, 8);
                   zeros(4), eye(4), zeros(4)];

if max(max(abs(CompanionMatrix-oo_.var.toto.CompanionMatrix)))>1e-12
   error('Companion matrix is wrong.')
end
