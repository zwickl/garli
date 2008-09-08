#!/bin/sh
nclv="ncl-2.1.04"
tar xfvz ${nclv}.tar.gz || exit
cd ${nclv} || exit
./configure --prefix=`pwd`/installed --disable-shared --enable-static --disable-dependency-tracking || exit
make || exit
make install || exit
make installcheck || exit
cd ..
./configure --with-ncl=`pwd`/${nclv}/installed || exit
make || exit
