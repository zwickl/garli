AM_CPPFLAGS = @INCLUDE@ @CPPFLAGS@
AM_CXXFLAGS = @CXXFLAGS@
AM_LDFLAGS  = @LDFLAGS@
 
bin_PROGRAMS = Garli@VERSION@

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
	hashdefines.h \
	individual.h \
	jph.h \
	linalg.h \
	memchk.h \
	mlhky.h \
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
	topologylist.h \
	translatetable.h \
	tree.h \
	treenode.h \
	utility.h 

Garli@VERSION@_SOURCES = \
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
	topologylist.cpp \
	translatetable.cpp \
	tree.cpp \
	treenode.cpp \
	mpitrick.cpp


Garli@VERSION@_LDADD =  $(LDADD) @GARLI_LIBS@
