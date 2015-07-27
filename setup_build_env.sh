function check()
{
	if test $? == 1
	then
		echo Failure: $aarch $arch $variant
		exit 1
	fi
}

echo XXX BUILD_AND_INSTALL XXX


variant=$1
threebdir=$2
TOON=$3
CVD=$4
GVARS=$5
threebversion=${threebdir}.tar.gz
J=-j8

	
echo $variant
echo $arch
pwd

CC=gcc	
AR=ar
LD=ld
RANLIB=ranlib
PIC=-fPIC
pref=
jni=

generate_arith=0

if [ x$variant == xmingw32 ]
then
	host="--host i686-w64-mingw32  --without-pthread"
	#Argh, the MingW *environment* provides the microsoft FPU control flags, but
	#The cross environment doesn't! Define them here...
	CC="i686-w64-mingw32-gcc"
	LAPACK_CCFLAGS="-DUSE_CLOCK -D_EM_DENORMAL=0x00080000 -D_EM_UNDERFLOW=0x00000002 -D_EM_INEXACT=0x00000001  -D_MCW_EM=0x0008001F"

	
	AR=i686-w64-mingw32-ar 
	prefdir=/tmp/mingw32
	CVD_LIBS="-L$prefdir/lib -lz"

	RANLIB=i686-w64-mingw32-ranlib
	pref="--prefix=/tmp/mingw32"
	jni=--with-jni=/tmp/java-win
	export CPPFLAGS=-I/tmp/mingw32/include
	export LDFLAGS=-L/tmp/mingw32/lib
	
	PIC=
	generate_arith=1
elif [ x$variant == xmingw64 ]
then
	host="--host x86_64-w64-mingw32  --without-pthread"
	#Argh, the MingW *environment* provides the microsoft FPU control flags, but
	#The cross environment doesn't! Define them here...
	CC="x86_64-w64-mingw32-gcc"
	LAPACK_CCFLAGS="-DUSE_CLOCK -D_EM_DENORMAL=0x00080000 -D_EM_UNDERFLOW=0x00000002 -D_EM_INEXACT=0x00000001  -D_MCW_EM=0x0008001F"

	prefdir=/tmp/mingw64
	CVD_LIBS="-L$prefdir/lib -lz"

	AR=x86_64-w64-mingw32-ar 
	RANLIB=x86_64-w64-mingw32-ranlib
	pref="--prefix=$prefdir"
	
	jni=--with-jni=/tmp/java-win
	export CPPFLAGS=-I/tmp/mingw64/include
	export LDFLAGS=-L/tmp/mingw64/lib
	
	PIC=
	generate_arith=1
fi


if [ $generate_arith == 1 ]
then
	#Machine dependent paramsfor CLAPACK
cat > /tmp/CLAPACK-3.2.1/F2CLIBS/libf2c/arith.h << FOO
#define IEEE_8087
#define Arith_Kind_ASL 1
#define Double_Align
#define QNaN0 0x0
#define QNaN1 0xfff80000
FOO

fi





