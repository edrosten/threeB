. /tmp/setup_build_env.sh


echo AAe Compiling 3B
cd /tmp/$threebdir
export CXXFLAGS="-DTOON_CLAPACK -O3"

#./configure $host --with-imagej=/tmp/ImageJ/ij.jar $jni
#make clean
#make $J
#touch ThreeBRunner.h


if [ x$variant != xstatic  ]
then
	./configure $host --with-imagej=/tmp/ImageJ/ij.jar $jni
	make clean
	#make $J
	#touch ThreeBRunner.h

	if [ x$variant == xmingw32 ] || [ x$variant == xmingw64 ]
	then
		make $J threeB_jni.dll
	else
		make $J libthreeB_jni.so
		make $J three_B.jar
		echo AAx result: $?
	fi
else
	./configure $host
	make clean
	make $J multispot5_static
fi

exit 0
