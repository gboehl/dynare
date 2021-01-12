[name='eq:x3']
x3 = a2*x3(-1)+b3*nu3+c3*y3(-2)+ex3;

[name='eq:y3']
diff(-diff(y3)) = d3*diff(x3(-1))+e3*nu3(-1)+f3*z3(-1)+ey3;

[name='eq:z3']
diff(log(-z3)) = g3*log(-z3(-1))+ez3;