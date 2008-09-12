#!/bin/sh
set -x
nl=`ls -l ncl*z | wc -l`
if ! test $nl = 1
then
	echo "You have more than one ncl version..."
	nclv=`ls ncl*gz | tail -n1 | sed -E 's/.tar.gz//'`
	echo "Using most recent:  $nclv"
else
	nclv=`ls ncl*z`
fi
if ! test -d ${nclv}
then
	tar xfvz ${nclv}.tar.gz || exit
fi
cd ${nclv} || exit
./configure --prefix=`pwd`/installed --disable-shared --enable-static || exit
make || exit
make install || exit
make installcheck || exit
cd ..
./configure $@ --prefix=`pwd` --with-ncl=`pwd`/${nclv}/installed || exit
make || exit
make install || exit
