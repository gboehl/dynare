function b=test_constraint(x,constr)
b = zeros(size(x,1),size(constr,1));
for i=1:size(constr,1)
    switch constr{i,4}
      case '<'
        b(:,i) = ~(x(:,constr{i,2}) < constr{i,3});
      case '<='
        b(:,i) = ~(x(:,constr{i,2}) <= constr{i,3});
      case '>'
        b(:,i) = ~(x(:,constr{i,2}) > constr{i,3});
      case '>='
        b(:,i) = ~(x(:,constr{i,2}) >= constr{i,3});
      otherwise
        error('OCCBIN: wrong inequality sign')
    end
end