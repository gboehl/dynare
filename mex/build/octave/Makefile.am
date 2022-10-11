ACLOCAL_AMFLAGS = -I ../../../m4

SUBDIRS = mjdgges kronecker bytecode block_kalman_filter sobol perfect_foresight_problem num_procs block_trust_region disclyap_fast libkordersim local_state_space_iterations folded_to_unfolded_dr k_order_simul k_order_mean cycle_reduction logarithmic_reduction

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

install-exec-local:
	$(MKDIR_P) $(DESTDIR)$(pkglibdir)/mex/octave
