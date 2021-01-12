function dy=extFunFirstDerivMatrix(a,b,c)
dy = [ c(1)^2 0        5*c(2)^2 0        2*a(1)*c(1) 10*b(1)*c(2);
       0      3*c(1)^2 0        7*c(2)^2 6*a(2)*c(1) 14*b(2)*c(2) ];
end
