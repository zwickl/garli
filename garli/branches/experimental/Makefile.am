ACLOCAL_AMFLAGS = -I config/m4

EXTRA_DIST = \
	build_garli.sh \
	README.txt \
	QuickStart.txt \
	example \
	project \
	doc \
	tests\
	ncl-2.1.18.tar.gz

SUBDIRS = src tests

dist-hook:
	find "$(distdir)/doc" -depth -name .svn -and -type d -and -exec rm -rf {} \; 
	find "$(distdir)/project" -depth -name .svn -and -type d -and -exec rm -rf {} \; 
	find "$(distdir)/example" -depth -name .svn -and -type d -and -exec rm -rf {} \; 
	find "$(distdir)/tests" -depth -name .svn -and -type d -and -exec rm -rf {} \; 
