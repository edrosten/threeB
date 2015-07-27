. /tmp/setup_build_env.sh

set -o xtrace

echo AAe Compiling 3B
cd /tmp/$threebdir
export CXXFLAGS="-DTOON_CLAPACK -O3"

#./configure $host --with-imagej=/tmp/ImageJ/ij.jar $jni
#make clean
#make $J
#touch ThreeBRunner.h


if [ x$variant != xstatic  ]
then
	#make $J
	#touch ThreeBRunner.h

	if [ x$variant == xmingw32 ] || [ x$variant == xmingw64 ]
	then
		LDFLAGS="$LDFLAGS -lpng -ltiff -ljpeg -lz" LIBS="-lpng -ltiff -ljpeg -lz" ./configure $host --with-imagej=/tmp/ImageJ/ij.jar $jni
		make clean
		make $J threeB_jni.dll multispot5_static
	else
		./configure $host --with-imagej=/tmp/ImageJ/ij.jar $jni
		make clean
		make $J libthreeB_jni.so
		make $J three_B.jar
		echo AAx result: $?
	fi
else
	./configure $host
	make clean
	make $J multispot5_static make_grid_markup_static

fi

exit 0
