
check-local:
	$(srcdir)/runtests.sh $(srcdir) $(top_builddir)/src/Garli$(EXEEXT) @NCL_BIN_DIR@/NEXUSvalidator
gpucheck:
	$(srcdir)/runtests.sh $(srcdir) $(top_builddir)/src/Garli$(EXEEXT) @NCL_BIN_DIR@/NEXUSvalidator -F GPU
dblgpucheck:
	$(srcdir)/runtests.sh $(srcdir) $(top_builddir)/src/Garli$(EXEEXT) @NCL_BIN_DIR@/NEXUSvalidator -F GPU DOUBLE
snglcpucheck:
	$(srcdir)/runtests.sh $(srcdir) $(top_builddir)/src/Garli$(EXEEXT) @NCL_BIN_DIR@/NEXUSvalidator -F CPU SINGLE
dblcpucheck:
	$(srcdir)/runtests.sh $(srcdir) $(top_builddir)/src/Garli$(EXEEXT) @NCL_BIN_DIR@/NEXUSvalidator -F CPU DOUBLE
ssecpucheck:
	$(srcdir)/runtests.sh $(srcdir) $(top_builddir)/src/Garli$(EXEEXT) @NCL_BIN_DIR@/NEXUSvalidator -F CPU SSE
