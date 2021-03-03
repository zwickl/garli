
if BEAGLE 
check-local:
	$(srcdir)/runtests.BEAGLE.sh	$(srcdir) $(top_builddir)/src/Garli$(EXEEXT) @NCL_BIN_DIR@/NEXUSvalidator
gpucheck:
	$(srcdir)/runtests.BEAGLE.sh $(srcdir) $(top_builddir)/src/Garli$(EXEEXT) @NCL_BIN_DIR@/NEXUSvalidator -F GPU
dblgpucheck:
	$(srcdir)/runtests.BEAGLE.sh $(srcdir) $(top_builddir)/src/Garli$(EXEEXT) @NCL_BIN_DIR@/NEXUSvalidator -F GPU DOUBLE
snglcpucheck:
	$(srcdir)/runtests.BEAGLE.sh $(srcdir) $(top_builddir)/src/Garli$(EXEEXT) @NCL_BIN_DIR@/NEXUSvalidator -F CPU SINGLE
dblcpucheck:
	$(srcdir)/runtests.BEAGLE.sh $(srcdir) $(top_builddir)/src/Garli$(EXEEXT) @NCL_BIN_DIR@/NEXUSvalidator -F CPU DOUBLE
ssecpucheck:
	$(srcdir)/runtests.BEAGLE.sh $(srcdir) $(top_builddir)/src/Garli$(EXEEXT) @NCL_BIN_DIR@/NEXUSvalidator -F CPU SSE
else
check-local:
	$(srcdir)/runtests.sh $(srcdir) $(top_builddir)/src/Garli$(EXEEXT) @NCL_BIN_DIR@/NEXUSvalidator
endif

