AM_CPPFLAGS = @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@
 
bin_PROGRAMS = Garli

EXTRA_DIST = threadfunc.cpp \
	mpifuncs.cpp 

noinst_HEADERS = \
	adaptation.h \
	bipartition.h \
	clamanager.h \
	condlike.h \
	configoptions.h \
	configreader.h \
	calculationmanager.h \
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

if BEAGLE 
Garli_SOURCES = \
	adaptation.cpp \
	bipartition.cpp \
	calculationmanager.cpp \
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
	optimization.cpp \
	population.cpp \
	rng.cpp \
	sequencedata.cpp \
	set.cpp \
	subsetcalculationmanager.cpp \
	translatetable.cpp \
	tree.cpp \
	treenode.cpp \
	mpitrick.cpp
else

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
	optimization.cpp \
	population.cpp \
	rng.cpp \
	sequencedata.cpp \
	set.cpp \
	translatetable.cpp \
	tree.cpp \
	treenode.cpp \
	mpitrick.cpp
endif

Garli_LDADD =  $(LDADD) @GARLI_LIBS@

install-exec-hook:
	cd $(DESTDIR)$(bindir) && \
	  mv -f Garli$(EXEEXT) Garli-$(VERSION)$(EXEEXT) && \
	  $(LN_S) Garli-$(VERSION)$(EXEEXT) Garli$(EXEEXT)

