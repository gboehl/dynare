MODFILES = \
	asset.mod \
	c20.mod \
	czech2.mod \
	dm7.mod \
	example1.mod \
	example1_optim.mod \
	gentay1a.mod \
	judd.mod \
	judd_norm.mod \
	kp1980_1.mod \
	kp1980_2.mod \
	kp1980_3.mod \
	lucas78.mod \
	m_0_3_0_0_0_0_0_0.mod \
	m_1_3_0_0_0_0_0_0.mod \
	m_1_3_0_0_0_0_0_1.mod \
	or0a.mod \
	portfolio4.mod \
	portfolio4_norm.mod \
	portfolio.mod \
	psd_exo3.mod \
	q3a2.mod \
	q3a50.mod \
	sedmodel1.mod \
	swma_pie.mod \
	test1.mod \
	test2a.mod \
	test2.mod \
	test3.mod \
	test4.mod \
	test5.mod \
	test6.mod \
	test7.mod \
	test.mod

EXTRA_DIST = $(MODFILES) \
	sw_euro.mod # This one crashes at steady state computation

check-local: $(MODFILES:%.mod=%.jnl)

%.jnl: %.mod
	../src/dynare++ --sim 2 $<

clean-local:
	rm -f *.jnl *_f.m *_ff.m *.dump
