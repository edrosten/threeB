base=/chroot
mirror=http://old-releases.ubuntu.com/ubuntu/
distro=precise
branch=`git branch | awk '/\*/{print $2}'`

if ! git status | awk '$1=="modified:"{exit 1}'
then
	echo "There are uncommitted changes. This means that the build"
	echo "Will not represent a sane version. Please commit all changes"
	echo "before building."
	exit 7
fi


doxyproj=`awk '$1=="PROJECT_NUMBER"{print $3}' Doxyfile`
autoproj=`awk -F'[,)]' '/AC_INIT/{print $2}' configure.ac`

if [ "$doxyproj" != "$autoproj" ]
then
	echo "Fatal Error: inconsistent project number!"
	echo "Doxygen: $doxyproj"
	echo "Autoconf: $autoproj"
	echo "Please update the project number in one or both of those files."
fi

git_hash=`git rev-parse HEAD`

if echo $autoproj | grep -q dev
then
	autoproj="$autoproj-`echo $git_hash| head -c 6`"
fi

threebdir=threeB-$autoproj
threebversion=${threebdir}.tar.gz
downloaddir=$PWD/build-downloads/


echo "Starting 3B build"
echo "Building verson $threebdir"

echo "Creating tar ball..."

#######################
#
# Now create a tar gz file. 
# Shove in a version specific About file.
#

mkdir "$threebdir" || exit 9
mkdir "$threebdir"/jar || exit 9

cp jar/about.txt "$threebdir/jar"
echo "Version: $threebdir" >> "$threebdir/jar/about.txt"
echo "git hash: $git_hash" >> "$threebdir/jar/about.txt"

echo "#define BUILDVERSION  \"$threebdir\"" >> "$threebdir/version.cc"
echo "#define BUILDHASH     \"$git_hash\"" >> "$threebdir/version.cc"
echo "#define BUILDDATE     \"$(date +'%Y:%m:%d-%T%z %s')\""  >> "$threebdir/version.cc"
echo "#define BUILDHOST     \"$(uname -a)\"" >> "$threebdir/version.cc"


git archive --prefix $threebdir/ $git_hash  > $threebdir.tar
tar -rf $threebdir.tar "$threebdir/jar/about.txt" 
tar -rf $threebdir.tar "$threebdir/version.cc"
gzip -9 $threebdir.tar

if ! [ -e $threebversion ]
then
	echo Error: $threebversion missing.
	exit 5
fi


CVDNAME=libcvd-20121025
GVARSNAME=gvars-3.0
TOONNAME=TooN-2.0

CVD=$CVDNAME.tar.gz
TOON=$TOONNAME.tar.gz
GVARS=$GVARSNAME.tar.gz



set -ex
mkdir -p $base

function check()
{
	if test $? != 0
	then
		echo Failure: $aarch $arch $variant
		exit
	fi
}
export arch
export aarch
export variant


#Fetch all the required source files, if they're not already here.
mkdir -p $downloaddir
(
	cd $downloaddir

	[ -f $TOON ] || wget http://www.edwardrosten.com/cvd/$TOON
	check

	[ -f $CVD ] || wget http://www.edwardrosten.com/cvd/$CVD
	check

	[ -f $GVARS ] || wget http://www.edwardrosten.com/cvd/$GVARS 
	check

	[ -f clapack.tgz ] || 	wget http://www.netlib.org/clapack/clapack.tgz
	check

	[ -f ij145.zip ] || wget http://rsbweb.nih.gov/ij/download/zips/ij145.zip
	check



	[ -f jpegsrc.v9a.tar.gz ] || wget http://www.ijg.org/files/jpegsrc.v9a.tar.gz
	check

	[ -f libpng-1.6.18.tar.xz  ] ||  wget https://download.sourceforge.net/libpng/libpng-1.6.18.tar.xz
	check

	[ -f tiff-3.9.7.tar.gz ] || wget https://download.osgeo.org/libtiff/tiff-3.9.7.tar.gz
	check

	[ -f zlib-1.2.8.tar.gz ] || wget http://prdownloads.sourceforge.net/libpng/zlib-1.2.8.tar.gz
	check


)


#list="i386 amd64 amd64_mingw32 amd64_mingw64 amd64_static i386_static"

#Create an output directory
dist=dist-$branch-$git_hash/ThreeB
mkdir -p $dist
check

function execute_build()
{
	#Order of things
	#
	# First create i386_base and amd64_base if they do not exist.
	# These are base systems.
	#
	# Then create ${list}_dev if they do not exist. Start by copying
	# the appropriate base system. These are pristine development
	# environments
	#
	# Then, copy ${list}_dev to ${list}_build
	# And build all but 3B. These are 3B build environments with all
	# but 3B built. Start by copying these from ${list}_dev
	#
	# Then nuke ${list}_3b, copy ${list}_built to ${list}_3b
	# This makes the _3b directory a full system with all the custom
	# libraries built. All that is required is to then copy over and
	# build 3B.
	#
	# The build of 3B is very well behaved and doesn't ever
	# escape its own directory, so we can probably avoid the rather slow
	# copy which now dominates the build time. I.e. ro rebuild, if it exists
	# just nuke the 3B directory inside _build.

	#First, create the base images if they do not exist.
	echo "Creating base systems..."
	arch=amd64  # i386
	target=$base/${distro}_${arch}_base


	if ! [ -e "$target" ]
	then
		echo "Need to install " $target
		sleep 1

		sudo debootstrap --no-check-gpg --variant=buildd --arch $arch ${distro} $target $mirror
		check

		sudo mount -o bind /proc $target/proc
		check

		sudo cp /etc/resolv.conf $target/etc/resolv.conf
		check

		sudo chroot $target locale-gen en_GB.UTF-8
		check

		sudo chroot $target apt-get install -y  --force-yes openjdk-6-jre-headless
		check 

		sudo chroot $target apt-get install -y  --force-yes openjdk-6-jdk cvs wget zip vim
		check

		sudo umount $target/proc
		check

	else
		echo "$target already exists. Yay!"	
	fi

	echo "Done with base system"
	sleep 1

	echo "Creating the individual development systems..."

	for aarch in $list
	do
		{
			arch=${aarch%_*}
			variant=${aarch#*_}
			target=$base/${distro}_${aarch}_dev
			source=$base/${distro}_${arch}_base

			if ! [ -e "$target" ]
			then
				
				echo YYY INSTALLING $aarch
				echo Getting $aarch from $source

				sleep 1

				cp -a "$source" "$target"

				#sudo mount -o bind /proc $target/proc
				#check



				#Build static exe, rather than plugin	
				if [ x$variant == xstatic ]
				then
					echo deb $mirror ${distro} main restricted universe | sudo tee $target/etc/apt/sources.list.d/asdf.list

					check
					sudo chroot $target apt-get update
					check
					sudo chroot $target apt-get install -y --force-yes libjpeg-dev libpng-dev libtiff-dev 
					check
				fi

				if [ x$variant == xmingw32 ]
				then
					echo deb $mirror ${distro} main restricted universe | sudo tee $target/etc/apt/sources.list.d/asdf.list
					check
					sudo chroot $target apt-get update
					check
					sudo chroot $target apt-get install -y --force-yes mingw-w64 g++-mingw-w64-i686
					check
				fi

				if [ x$variant == xmingw64 ]
				then
					echo deb $mirror ${distro} main restricted universe | sudo tee $target/etc/apt/sources.list.d/asdf.list
					check
					sudo chroot $target apt-get update
					check
					sudo chroot $target apt-get install -y --force-yes mingw-w64 g++-mingw-w64-x86-64
					check
				fi

				#sudo umount $target/proc
				#check

			else
				echo $target already exists. Yay.
			fi
		}& 
	done
	wait

	echo "Done creating individual development systems"
	sleep 1


	for aarch in $list
	do
		arch=${aarch%_*}
		variant=${aarch#*_}
		target=$base/${distro}_${aarch}_build
		source=$base/${distro}_${aarch}_dev

		if ! [ -e "$target" ]
		then
		
			echo Copying to $target
			sudo cp -a "${source}" "${target}"

			echo ZZZ GETTING $aarch $target
			sleep 1
			(
				cp build_and_install.sh \
				   setup_build_env.sh \
				   build_3b.sh \
				   clapack_mingw.patch \
				   $target/tmp
				check

				cd $target/tmp
				check
			
				if true
				then
					cp $downloaddir/$TOON .
					check

					cp $downloaddir/$CVD .
					check

					cp $downloaddir/$GVARS  .
					check

					cp $downloaddir/clapack.tgz .
					check

					cp $downloaddir/ij145.zip .
					check

					tar -xzf $TOON
					check

					tar -xzf $GVARS 
					check

					tar -xzf $CVD
					check
				
					tar -xzf clapack.tgz
					check

					unzip ij145.zip


					if [ x$variant == xmingw32 ] || [ x$variant == xmingw64 ]
					then
						cp $downloaddir/jpegsrc.v9a.tar.gz .
						check
						tar -xzf jpegsrc.v9a.tar.gz 
						check
						
						cp $downloaddir/libpng-1.6.18.tar.xz .
						check
						tar -xJf libpng-1.6.18.tar.xz
						check

						cp $downloaddir/tiff-3.9.7.tar.gz .
						check
						tar -xzf tiff-3.9.7.tar.gz
						check
						
						cp $downloaddir/zlib-1.2.8.tar.gz .
						check
						tar -xvf zlib-1.2.8.tar.gz
						check

						mkdir java-win
						check
						cp ../usr/lib/jvm/java-6-openjdk-amd64/include/jni.h java-win
						check
						
						#The JNI interface is very, VERY stable. Here is the machine dependent
						#file for 32 bit java 1.7.0.2 on Win32
						#
						# Note that the Java types are the same length on any platform,
						# i.e. jlong is _always_ 64 bit.
						#
						# The Windows API is also quite stable, i.e. long is
						# 32 bits on 32 or 64 bit systems
						cat > java-win/jni_md.h << FOO

#ifndef _JAVASOFT_JNI_MD_H_
#define _JAVASOFT_JNI_MD_H_

#define JNIEXPORT __declspec(dllexport)
#define JNIIMPORT __declspec(dllimport)
#define JNICALL __stdcall

typedef long jint;
typedef __int64 jlong;
typedef signed char jbyte;

#endif /* !_JAVASOFT_JNI_MD_H_ */
FOO
					fi
				fi


				sudo chroot $target bash /tmp/build_and_install.sh $variant $threebdir $TOONNAME $CVDNAME $GVARSNAME
				check
			)
		else
			echo "Build env $target already exists."
		fi
	done


	for aarch in $list
	do
		arch=${aarch%_*}
		variant=${aarch#*_}
		target=$base/${distro}_${aarch}_build
		
		echo Nuking 3B in $target
		sudo rm -rf "$target"/tmp/$threebdir

		#This is astonishingly required to use java under ubuntu, now!
		sudo mount -o bind /proc $target/proc
		check
		(

			cp $threebversion $target/tmp
			check

			cp build_and_install.sh \
			   setup_build_env.sh \
			   build_3b.sh \
			   clapack_mingw.patch \
			   $target/tmp
			check

			cd $target/tmp
			check
			tar -xzf $threebversion
			check
			
			sudo chroot $target bash /tmp/build_3b.sh $variant $threebdir $TOONNAME $CVDNAME $GVARSNAME
			check
			echo done1
		)

		check
		sudo umount $target/proc
		echo done2

	done



	if echo  $list | tr ' ' '\n' | grep -q amd64_static
	then
		cp $base/${distro}_amd64_static_build/tmp/$threebdir/multispot5_static $dist/multispot5_static_amd64_linux
		cp $base/${distro}_amd64_static_build/tmp/$threebdir/make_grid_markup_static $dist/make_grid_markup_static_amd64_linux
		check
	fi

	if echo  $list | tr ' ' '\n' | grep -q i386_static
	then
		cp $base/${distro}_i386_static_build/tmp/$threebdir/multispot5_static $dist/multispot5_static_i386_linux
		cp $base/${distro}_i386_static_build/tmp/$threebdir/make_grid_markup_static $dist/make_grid_markup_static_i386_linux
		check
	fi

	if echo  $list | tr ' ' '\n' | grep -q 'i386$'
	then
		cp $base/${distro}_i386_build/tmp/$threebdir/libthreeB_jni.so $dist/libthreeB_jni_32.so
		check
		cp $base/${distro}_i386_build/tmp/$threebdir/three_B.jar $dist/
		check
	fi

	if echo  $list | tr ' ' '\n' | grep -q 'amd64$'
	then
		cp $base/${distro}_amd64_build/tmp/$threebdir/libthreeB_jni.so $dist/libthreeB_jni_64.so
		check
		cp $base/${distro}_amd64_build/tmp/$threebdir/three_B.jar $dist/
		check
	fi


	if echo  $list | tr ' ' '\n' | grep -q 'mingw32$'
	then
		cp $base/${distro}_amd64_mingw32_build/tmp/$threebdir/threeB_jni.dll $dist/threeB_jni_32.dll
		cp $base/${distro}_amd64_mingw32_build/tmp/$threebdir/multispot5_static $dist/multispot5_static_i386_win32.exe
		check
	fi

	if echo  $list | tr ' ' '\n' | grep -q 'mingw64$'
	then
		cp $base/${distro}_amd64_mingw64_build/tmp/$threebdir/threeB_jni.dll $dist/threeB_jni_64.dll
		cp $base/${distro}_amd64_mingw64_build/tmp/$threebdir/multispot5_static $dist/multispot5_static_amd64_win64.exe
		check
	fi
}

export distro=precise
#export list="i386 amd64_mingw64" #i386 required for JAR file
export list="amd64 amd64_mingw64"
execute_build

export list="amd64_static"
export distro=lucid
execute_build


#export distro=precise
#export list="amd64_mingw64"
#execute_build


echo $dist
cd $dist/..
zip -r ThreeB.zip ThreeB/
tar -cvzf  ThreeB.tar.gz ThreeB/
tar -cvjf  ThreeB.tar.bz2 ThreeB/
tar -cvJf  ThreeB.tar.xz ThreeB/
cd ..


echo $dist
