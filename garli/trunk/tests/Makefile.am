
if BEAGLE 
check-local:
	$(srcdir)/runtests.BEAGLE.sh	$(srcdir) $(top_builddir)/src/Garli$(EXEEXT) @NCL_BIN_DIR@/NEXUSvalidator
else
check-local:
	$(srcdir)/runtests.sh $(srcdir) $(top_builddir)/src/Garli$(EXEEXT) @NCL_BIN_DIR@/NEXUSvalidator
endif

