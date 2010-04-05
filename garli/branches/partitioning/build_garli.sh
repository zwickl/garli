#!/bin/sh

if [ $# -gt 0 ]
then
        if [ "$1" = "--ncl-svn" ]
        then
		shift # this shifts the first cmd line argument out so that the rest can be passed to GARLI configure
                echo "***CHECKING OUT NCL LIBRARY SOURCE VIA SUBVERSION***"
                svn co http://ncl.svn.sourceforge.net/svnroot/ncl/branches/v2.1 ncl-svn || exit
                nclv="ncl-svn"
                cd ${nclv} || exit
                sh bootstrap.sh || exit
                cd ..
        elif [ "$1" = "-h" ] || [ "$1" = "--help" ]
	then
                echo "Usage ./$0 [--svn-ncl] [-h] [arguments to GALRI configure script]"
                echo "     --ncl-svn   Check out current NCL source via anonymous svn, build NCL, then GARLI"
                echo "     --ncl-dist  Automatically build NCL from a ncl-2.1.xx.tar.gz distribution"
                echo "                 in this directory, then build GARLI (default)"
                echo "  -h --help      Output this help and exit"
                echo "  [other args]   Other arguments are passed to GARLI's configure invocation"
                echo
                exit
        fi
fi

#if NCL wasn't checked out above
if [ -z "${nclv}" ]
then
	if [ "$1" = "--ncl-dist" ]
	then
		shift # this shifts the first cmd line argument out so that the rest can be passed to GARLI configure
	fi
	echo "***BUILDING NCL LIBRARY FROM SOURCE DISTRIBUTION***"
	nl=`ls -l ncl*.gz | wc -l`
	if [ $nl -eq 0 ]
	then
		echo "ERROR: No ncl-2.1.xx.tar.gz distributions found."
		echo "  Provide one or try \"$0 --ncl-svn\" to checkout NCL via subversion and automatically build NCL and GARLI."
		exit
	elif [ ! $nl -eq 1 ]
	then
		echo "You have more than one NCL version..."
        	nclv=`ls ncl*.gz | tail -n1 | sed 's/.tar.gz//'`
        	echo "Using most recent:  $nclv"
	else
        	nclv=`ls ncl*.gz | sed 's/.tar.gz//'`
	fi
	if [ ! -d ${nclv} ]
	then
       		tar xfvz ${nclv}.tar.gz || exit
	fi
fi
cd ${nclv} || exit
echo "CONFIGURING NCL ..."
./configure --prefix=`pwd`/installed --disable-shared --enable-static || exit
make || exit
echo "BUILDING NCL ..."
make install || exit
make installcheck || exit
cd ..

echo "CONFIGURING GARLI ..."
if [ ! -f configure ]
then
	if [ -f bootstrap.sh ]
	then
		sh bootstrap.sh || exit
	else
		echo "Neither configure nor bootstrap.sh found.  This is not a complete distribution."
	fi
fi
./configure $@ --prefix=`pwd` --with-ncl=`pwd`/${nclv}/installed || exit
echo "BUILDING GARLI ..."
make || exit
make install || exit
