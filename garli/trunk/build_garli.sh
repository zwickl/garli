#!/bin/sh


if [ $# -gt 0 ]
then
        if [ "$1" = "--ncl-svn" -o "$1" = "--ncl-sourceforge" ]
        #if [ "$1" = "--ncl-svn" ]
        then
            if [ -d ncl-svn ]
            then
                echo "***NCL LIBRARY SOURCE FROM SUBVERSION ALREADY EXISTS***"
                echo "***CURRENT COPY WILL BE USED AS-IS.  UPDATE IT MANUALLY OR***"
                echo "***DELETE THE ncl-svn DIRECTORY TO GET THE LATEST NCL SOURCE***"
            else
                echo "***CHECKING OUT NCL LIBRARY SOURCE VIA SUBVERSION***"
                if [ "$1" = "--ncl-svn" ];then
                    svn co https://github.com/mtholder/ncl/trunk ncl-svn || exit
                else
                    svn co http://svn.code.sf.net/p/ncl/code/branches/v2.1 ncl-svn || exit
                fi
            fi
            shift # this shifts the first cmd line argument out so that the rest can be passed to GARLI configure
            nclv="ncl-svn"
            cd ${nclv} || exit
            sh bootstrap.sh || exit
            cd ..
        elif [ "$1" = "--ncl-git" ]
        then    
            TARG="ncl-git"
            nclv=$TARG
            if [ -d $TARG ] 
            then
                echo "***NCL LIBRARY SOURCE FROM GIT ALREADY EXISTS***"
                echo "***CURRENT COPY WILL BE USED AS-IS.  UPDATE IT MANUALLY OR***"
                echo "***DELETE THE $TARG  DIRECTORY TO GET THE LATEST NCL SOURCE***"
            else
                echo "***CHECKING OUT NCL LIBRARY SOURCE VIA GIT***"
                    git clone https://github.com/mtholder/ncl.git $TARG|| exit
            fi
            shift # this shifts the first cmd line argument out so that the rest can be passed to GARLI configure
            cd ${TARG} || exit
            sh bootstrap.sh || exit
            cd ..
        elif [ "$1" = "-h" ] || [ "$1" = "--help" ]
        then
                echo "Usage ./$0 [--svn-ncl] [-h] [arguments to GALRI configure script]"
                echo "     --ncl-svn            Check out current NCL v2.1 source from github repo via anonymous svn, build NCL, then GARLI"
                echo "                              (prefer over ncl-sourceforge)"
                echo "     --ncl-sourceforge    Check out current NCL v2.1 source from sourceforge repo via anonymous svn,"
                echo "                              build NCL, then GARLI"
                echo "     --ncl-dist           Automatically build NCL from a ncl-2.1.xx.tar.gz distribution"
                echo "                              in this directory, then build GARLI (default)"
                echo "  -h --help               Output this help and exit"
                echo "  [other args]            Other arguments are passed to GARLI's configure invocation"
                echo
                exit
        fi
fi

#if NCL wasn't checked out above
if [ -z "${TARG}" ]
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
env CXXFLAGS=-DNCL_CONST_FUNCS ./configure --prefix=`pwd`/installed --disable-shared --enable-static || exit
make || exit
echo "BUILDING NCL ..."
make install || exit
#make installcheck || exit
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
cp ${nclv}/example/gapcode/NEXUSgapcode bin/

