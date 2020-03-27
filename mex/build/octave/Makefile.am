ACLOCAL_AMFLAGS = -I ../../../m4

SUBDIRS = mjdgges kronecker bytecode block_kalman_filter sobol perfect_foresight_problem num_procs block_trust_region

# libdynare++ must come before gensylv, k_order_perturbation, dynare_simul_
if ENABLE_MEX_DYNAREPLUSPLUS
SUBDIRS += libdynare++ gensylv k_order_perturbation dynare_simul_ local_state_space_iterations
endif

if ENABLE_MEX_MS_SBVAR
SUBDIRS += ms_sbvar
endif

if ENABLE_MEX_KALMAN_STEADY_STATE
SUBDIRS += kalman_steady_state
endif

install-exec-local:
	$(MKDIR_P) $(DESTDIR)$(pkglibdir)/mex/octave
