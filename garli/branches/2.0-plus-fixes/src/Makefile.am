AM_CPPFLAGS = @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@
 
bin_PROGRAMS = Garli

EXTRA_DIST = mpitrick.cpp \
	threadfunc.cpp

noinst_HEADERS = \
	adaptation.h \
	bipartition.h \
	clamanager.h \
	condlike.h \
	configoptions.h \
	configreader.h \
	datamatr.h \
	defs.h \
	errorexception.h \
	funcs.h \
	garlireader.h \
	individual.h \
	linalg.h \
	memchk.h \
	model.h \
	mpifuncs.h \
	optimizationinfo.h \
	outputman.h \
	population.h \
	reconnode.h \
	rng.h \
	sequencedata.h \
	set.h \
	stopwatch.h \
	threaddcls.h \
	translatetable.h \
	tree.h \
	treenode.h \
	utility.h 

Garli_SOURCES = \
	adaptation.cpp \
	bipartition.cpp \
	condlike.cpp \
	configoptions.cpp \
	configreader.cpp \
	datamatr.cpp \
	funcs.cpp \
	garlimain.cpp \
	garlireader.cpp \
	individual.cpp \
	linalg.cpp \
	model.cpp \
	mpifuncs.cpp \
	optimization.cpp \
	population.cpp \
	rng.cpp \
	sequencedata.cpp \
	set.cpp \
	translatetable.cpp \
	tree.cpp \
	treenode.cpp \
	mpitrick.cpp


Garli_LDADD =  $(LDADD) @GARLI_LIBS@

install-exec-hook:
	cd $(DESTDIR)$(bindir) && \
	  mv -f Garli$(EXEEXT) Garli-$(VERSION)$(EXEEXT) && \
	  $(LN_S) Garli-$(VERSION)$(EXEEXT) Garli$(EXEEXT)

