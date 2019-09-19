OCTAVE ?= octave-cli
MATLAB ?= $(shell which matlab)

all: check-octave check-matlab

check-octave:
	@cd test ;\
	$(OCTAVE) --no-init-file --silent --no-history --eval "addpath([pwd() '/../dseries/src']); addpath([pwd() '/../src']); runtest;" && [ -f my/report/dir/report.pdf ]

check-matlab:
	@$(MATLAB) -nosplash -nodisplay -r "addpath([pwd '/dseries/src']); addpath([pwd '/src']); cd test; runtest; quit" && [ -f test/my/report/dir/report.pdf ]
