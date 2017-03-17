OCTAVE=octave-cli
MATLAB=`which matlab`

all: check-octave check-matlab

check-octave:
	@cd test ;\
	$(OCTAVE) --silent --no-history runtest.m

check-matlab:
	@$(MATLAB) -nosplash -nodisplay -r "cd test; runtest; quit"
