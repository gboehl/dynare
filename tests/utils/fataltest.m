function fataltest(a,b,n)
  if max(max(abs(a-b))) > 1e-5
    error(['Test error in test ' int2str(n)])
  end
