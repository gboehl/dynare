function [X,w]=stroud_judd_7.5.8(d)
    
    E = eye(d);
    X = cell(2*d,1);
    m = 1;
    for i=1:d
        X{m} = E(:,i);
        m = m+1;
        X{m} = -E(:,i);
