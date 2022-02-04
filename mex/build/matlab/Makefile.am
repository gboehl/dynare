ACLOCAL_AMFLAGS = -I ../../../m4

SUBDIRS = mjdgges kronecker bytecode block_kalman_filter sobol perfect_foresight_problem num_procs block_trust_region disclyap_fast libkordersim local_state_space_iterations folded_to_unfolded_dr k_order_simul k_order_mean

# libdynare++ must come before gensylv and k_order_perturbation
if ENABLE_MEX_DYNAREPLUSPLUS
SUBDIRS += libdynare++ gensylv libkorder k_order_perturbation k_order_welfare
endif

if ENABLE_MEX_MS_SBVAR
SUBDIRS += ms_sbvar
endif

if ENABLE_MEX_KALMAN_STEADY_STATE
SUBDIRS += kalman_steady_state
endif

if HAVE_M2HTML
html-local:
	echo "addpath ${M2HTML}" > run_m2html.m
	echo "cd ../../.." >> run_m2html.m
	echo "m2html('mFiles', 'matlab', 'htmlDir', 'doc/m2html', 'recursive', 'on', 'globalHypertextLinks', 'on')" >> run_m2html.m
	echo "quit" >> run_m2html.m
	$(MATLAB)/bin/matlab -r run_m2html # Not using the -nodisplay option, because with m2html it confuses the terminal (need to do a "reset")

clean-local:
	rm -f run_m2html.m
	rm -rf ../../../doc/m2html
endif

install-exec-local:
	$(MKDIR_P) $(DESTDIR)$(pkglibdir)/mex/matlab

EXTRA_DIST = mex.def mexFunction-MacOSX.map
