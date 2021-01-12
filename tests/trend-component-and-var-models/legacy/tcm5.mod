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
           A341 A342 A343 A344
           B11 B12 B21 B22;

A1 = randn(4);
A2 = randn(4);
A3 = randn(4);

A1(3:4,:) = 0;
A2(3:4,:) = 0;
A3(3:4,:) = 0;

A1(1:2,3:4) = 0;
A2(1:2,3:4) = 0;
A3(1:2,3:4) = 0;

A1(3:4,3:4) = eye(2); // y3 and y4 are pure random walks.


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

B = rand(2);
B(1,1) = -B(1,1);
B(2,2) = -B(2,2);

B11 = B(1,1);
B12 = B(1,2);
B21 = B(2,1);
B22 = B(2,2);

trend_component_model(model_name=toto, eqtags=['eq:y4', 'eq:y1', 'eq:y2', 'eq:y3'], targets=['eq:y3', 'eq:y4']);

model;

[name='eq:y1']
diff(y1) = B11*(y1(-1)-y3(-1)) + B12*(y2(-1)-y4(-1)) + 
     A111*diff(y1(-1)) + A112*diff(y2(-1)) +
     A211*diff(y1(-2)) + A212*diff(y2(-2)) +
     A311*diff(y1(-3)) + A312*diff(y2(-3)) + e1; 


[name='eq:y4']
y4 = y4(-1) + e4;


[name='eq:y2']
diff(y2) = B21*(y1(-1)-y3(-1)) + B22*(y2(-1)-y4(-1)) +
     A121*diff(y1(-1)) + A122*diff(y2(-1)) +
     A221*diff(y1(-2)) + A222*diff(y2(-2)) +
     A321*diff(y1(-3)) + A322*diff(y2(-3)) + e2; 

[name='eq:y3']
y3 = y3(-1) + e3; 

end;

[EC, AR, T] = get_companion_matrix_legacy('toto');

if max(max(abs(EC-B)))>1e-12
   error('Error component matrix is wrong.')
end

if max(max(abs(AR(:,:,1)-A1(1:2,1:2))))>1e-12
   error('First order autoregressive matrix is wrong.')
end

if max(max(abs(AR(:,:,2)-A2(1:2,1:2))))>1e-12
   error('Second order autoregressive matrix is wrong.')
end

if max(max(abs(AR(:,:,3)-A3(1:2,1:2))))>1e-12
   error('Third order autoregressive matrix is wrong.')
end