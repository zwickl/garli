#!/bin/sh
set -x
nl=`ls -l ncl*.gz | wc -l`
if ! test $nl = 1
then
	echo "You have more than one ncl version..."
	nclv=`ls ncl*.gz | tail -n1 | sed 's/.tar.gz//'`
	echo "Using most recent:  $nclv"
else
	nclv=`ls ncl*.gz | sed 's/.tar.gz//'`
fi
if ! test -d ${nclv}
then
	tar xfvz ${nclv}.tar.gz || exit
fi
cd ${nclv} || exit
env CXXFLAGS=-DNCL_CONST_FUNCS ./configure --prefix=`pwd`/installed --disable-shared --enable-static || exit
make || exit
make install || exit
make installcheck || exit
cd ..
env CXXFLAGS=-DNCL_CONST_FUNCS ./configure $@ --prefix=`pwd` --with-ncl=`pwd`/${nclv}/installed || exit
make || exit
make install || exit
