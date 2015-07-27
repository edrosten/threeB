. /tmp/setup_build_env.sh



if [ x$variant == xmingw32 ] || [ x$variant == xmingw64 ]
then
	echo AA1 Building jpeglib
	cd /tmp/jpeg-9a
	check
	./configure $host $pref && make $J && make install 
	check

	echo AA2 Building zlib
	cd /tmp/zlib-1.2.8
	check


	#Zlib is speshul. Too "cool" to use the "crap" autoconf or something, so naturally
	#it's a pain to compile.
	RANLIB="$RANLIB" AR="$AR" CC="$CC" ./configure --static $pref && make libz.a $J && make install
	check

	echo AA3 Building png
	cd /tmp/libpng-1.6.18
	check
	LDFLAGS="-L$prefdir/lib -lz " ./configure $host $pref 
	make $J
	check
	make install 
	check

	echo AA4 Building libtiff
	cd /tmp/tiff-3.9.7
	check
	LDFLAGS="-L$prefdir/lib -lz " ./configure $host $pref
	check
	make $J
	check
	make install 
	check
fi




cd /tmp/CLAPACK-3.2.1
cp make.inc.example make.inc
patch make.inc <<FOO 
24c24
< CC        = gcc
---
> CC        = $CC $LAPACK_CCFLAGS
27c27
< CFLAGS    = -O3 -I\$(TOPDIR)/INCLUDE
---
> CFLAGS    = -O3 -I\$(TOPDIR)/INCLUDE $PIC
30c30
< NOOPT     = -O0 -I\$(TOPDIR)/INCLUDE
---
> NOOPT     = -O0 -I\$(TOPDIR)/INCLUDE $PIC
53c53
< ARCH     = ar
---
> ARCH     = $AR
55c55
< RANLIB   = ranlib
---
> RANLIB   = $RANLIB
FOO
check

##No need to remove main from libf2c because since we are linking a library
##apparently doesn't pull in main(). Hmm
#patch F2CLIBS/libf2c/Makefile < ../clapack.patch
#check

#Remove nasty things which clock cross compiling
patch -p1 < ../clapack_mingw.patch

make $J blaslib
check
make $J f2clib
check
cd INSTALL
make ilaver.o slamch.o dlamch.o lsame.o
check
echo > second.c
cd ..
make $J lapacklib
echo AAg results: $?
check

set -o xtrace
echo AAh Copying libraries	
if [ x$variant == xmingw32 ] || [ x$variant == xmingw64 ]
then
	mkdir -p $prefdir/lib
	cp blas_LINUX.a     $prefdir/lib/libblas.a
	check
	cp lapack_LINUX.a   $prefdir/lib/liblapack.a
	check
	cp F2CLIBS/libf2c.a $prefdir/lib/libf2c.a
else
	mkdir -p /usr/local/lib
	cp blas_LINUX.a /usr/local/lib/libblas.a
	check
	cp lapack_LINUX.a /usr/local/lib/liblapack.a
	check
	cp F2CLIBS/libf2c.a /usr/local/lib/libf2c.a
	check
fi


	
echo AAa Building and installing TooN2
cd /tmp/$TOON
check
./configure $host $pref && make $J && make install 
check

echo AAb Building and installing GVars3
cd /tmp/$GVARS
check
./configure $host $pref --without-head --without-lang && make $J && make install 
check


echo AAc Building and installing libcvd
cd /tmp/$CVD
check
LIBS="$CVD_LIBS" ./configure $host $pref --disable-fast7 --disable-fast8 --disable-fast9 --disable-fast10 --disable-fast11 --disable-fast12 && make $J && make install 
check

/sbin/ldconfig


echo AAd Perhaps doung some MinGW hacks
cd /tmp/mingw32/lib
for i in *.a
do 
	$RANLIB $i
done
cd /tmp/mingw64/lib
for i in *.a
do 
	$RANLIB $i
done

pwd

exit 0
