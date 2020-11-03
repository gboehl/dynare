// test ifndef to prevent regressions like #1747
@#define marco=1
@#ifndef marco
    def=0;
@#else
    def=1;
@#endif

if ~def
    error('ifndef not reached')
end

var x;

model;
x;
end;