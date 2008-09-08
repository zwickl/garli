ACLOCAL_AMFLAGS = -I config/m4

EXTRA_DIST = \
	build_garli.sh \
	README.txt \
	QuickStart.txt \
	example \
	project \
	doc \
	ncl-2.1.04.tar.gz

SUBDIRS = src 

dist-hook:
	find "$(distdir)/doc" -depth -and -name .svn -and -type d -and -exec rm -rf {} \; 
	find "$(distdir)/project" -depth -and -name .svn -and -type d -and -exec rm -rf {} \; 
	find "$(distdir)/example" -depth -and -name .svn -and -type d -and -exec rm -rf {} \; 


