var y1 y2 y3 y4;

varexo e1 e2 e3 e4;

parameters A111 A112
           A121 A122
           A211 A212
           A221 A222
           A311 A312
           A321 A322
           B11 B12 B21 B22;

A1 = randn(2);
A2 = randn(2);
A3 = randn(2);

A111 = A1(1,1);
A112 = A1(1,2);
A121 = A1(2,1);
A122 = A1(2,2);

A211 = A2(1,1);
A212 = A2(1,2);
A221 = A2(2,1);
A222 = A2(2,2);

A311 = A3(1,1);
A312 = A3(1,2);
A321 = A3(2,1);
A322 = A3(2,2);

B = rand(2);
B(1,1) = -B(1,1);
B(2,2) = -B(2,2);

B11 = B(1,1);
B12 = B(1,2);
B21 = B(2,1);
B22 = B(2,2);

trend_component_model(model_name=toto, eqtags=['eq:y1', 'eq:y2', 'eq:y3', 'eq:y4'], trends=['eq:y3', 'eq:y4']);

model;

[name='eq:y1']
diff(y1) = B11*(y1(-1)-y3(-1)) + B12*(y2(-1)-y4(-1)) + 
     A111*diff(y1(-1)) + A112*diff(y2(-1)) +
     A211*diff(y1(-2)) + A212*diff(y2(-2)) +
     A311*diff(y1(-3)) + A312*diff(y2(-3)) + e1; 


[name='eq:y2']
diff(y2) = B21*(y1(-1)-y3(-1)) + B22*(y2(-1)-y4(-1)) +
     A121*diff(y1(-1)) + A122*diff(y2(-1)) +
     A221*diff(y1(-2)) + A222*diff(y2(-2)) +
     A321*diff(y1(-3)) + A322*diff(y2(-3)) + e2; 

[name='eq:y3']
y3 = y3(-1) + e3; 

[name='eq:y4']
y4 = y4(-1) + e4;

end;

[EC, AR, T] = get_companion_matrix('toto');

if max(max(abs(EC-B)))>1e-12
   error('Error component matrix is wrong.')
end

A1fake = AR(:,:,1);

if max(max(abs(A1fake-A1)))>1e-12
   error('First order autoregressive matrix is wrong.')
end

if max(max(abs(AR(:,:,2)-A2)))>1e-12
   error('Second order autoregressive matrix is wrong.')
end

if max(max(abs(AR(:,:,3)-A3)))>1e-12
   error('Third order autoregressive matrix is wrong.')
end
