AC_INIT(threeB,1.2-dev)
AC_COPYRIGHT(Copyright Edward Rosten and Susan Cox (2011, 2012, 2013))


orig_CXXFLAGS="$CXXFLAGS"

AC_PROG_AWK
AC_PROG_CXX

AC_LANG(C++)

################################################################################
#
# Useful macros
#

dnl APPEND(var, value)
dnl This appends vale to a shell variable var
define(APPEND, [$1="$$1 $2"])

dnl CHECK_FOR_OPTION(option, code, [success], [fail])
dnl option corresponds to the first parameter for the corresponding
dnl AC_ARG_WITH. If the library is present, then have_option and CVD_HAVE_OPTION
dnl are set. If with_option is "", it will be set to yes.
dnl
dnl The test code shall do nothing on success, and set the shell variable a to 1
dnl on failure. For convenience, LIBS  and CXXFLAGS are checkpointed 
dnl automatically. On failure, they are restored automatically.
dnl
dnl On success, "option" is appended to options, and the code "success" is run
dnl On failure including --without) fail is run.
define(CHECK_FOR_OPTION,[
        if test "$have_$1" != yes
        then
                if test "$with_$1" == ""
                then
                        with_$1=yes
                fi

                save_LIBS="$LIBS"
                save_CXXFLAGS="$CXXFLAGS"
                a=
                $2
                
                if test x$a == x
                then
                        AC_DEFINE(CVD_HAVE_[]m4_toupper($1), 1)
                        AC_SUBST(have_$1, yes)
                        APPEND(options, $1)
                        $3
                else
                        LIBS="$save_LIBS"
                        CXXFLAGS="$save_CXXFLAGS"
                fi
        fi
        
        if test "$have_$1" != yes
        then
                :
                $4
        fi
        ])



dnl TEST_AND_SET_CXXFLAG(flag, [program], [run])
dnl
dnl This attempts to compile a and run program with a certain compiler flag.
dnl If no program is given, then the minimal C++ program is compiled, and 
dnl this tests just the validity of the compiler flag. 
dnl
define([TEST_AND_SET_CXXFLAG],[
	if test "$3" == ""
	then
		AC_MSG_CHECKING([if compiler flag $1 works])	
	else
		AC_MSG_CHECKING([$3])	
	fi
	save_CXXFLAGS="$CXXFLAGS"
	APPEND(CXXFLAGS, [$1])

	m4_if([$2],[],[pushdef(prog, [int main(){}])], [pushdef(prog, [$2])])

	m4_if([$4],[run],
		[AC_RUN_IFELSE([AC_LANG_SOURCE([prog])], [cvd_conf_test=1],[cvd_conf_test=0], [cvd_conf_test=0])],
		[AC_COMPILE_IFELSE([AC_LANG_SOURCE([prog])], [cvd_conf_test=1],[cvd_conf_test=0])]
	)
	
	
	popdef([prog])

	if test $cvd_conf_test = 1
	then
		AC_MSG_RESULT(yes)
		ts_success=yes
	else
		AC_MSG_RESULT(no)
		CXXFLAGS="$save_CXXFLAGS"
		ts_success=no
	fi
])


TEST_AND_SET_CXXFLAG(-Wall)
TEST_AND_SET_CXXFLAG(-Wextra)
TEST_AND_SET_CXXFLAG(-W)
TEST_AND_SET_CXXFLAG(-fPIC)

if ! echo $host | grep mingw
then
        TEST_AND_SET_CXXFLAG(-fPIC)
        if test $ts_success == no
        then
                AC_MSG_WARN(["I don't know how to compile position independent code. The Java plugin won't work"])
        fi
fi



################################################################################
#
# Check for TooN
#
toon=yes
AC_CHECK_HEADERS(TooN/TooN.h, [], [toon=no])

comp=no
AC_MSG_CHECKING([TooN])
AC_COMPILE_IFELSE([AC_LANG_SOURCE([
	#include <TooN/TooN.h>
	int main(){TooN::Vector<2> v;return 0;}
])], [comp=yes], [toon=no])

AC_MSG_RESULT($comp)

if test $toon = no
then
	AC_MSG_ERROR([A working copy of TooN is required for this program.
Get TooN here: http://www.edwardrosten.com/cvd/toon.html
])
fi


################################################################################
#
# Check for TaG
#
#tag=yes
#AC_CHECK_HEADERS(tag/printf.h, [], [tag=no])
#AC_MSG_CHECKING([tag])
#AC_COMPILE_IFELSE([
#	#include <tag/printf.h>
#	int main(){Printf("Hello, world.\n");}
#], [comp=yes], [tag=no])
#
#AC_MSG_RESULT($comp)
#
#if test $tag = no
#then
#	AC_MSG_ERROR([A working copy of tag is required for this program.
#Get tag here: http://www.edwardrosten.com/cvd/tag.html
#])
#fi
#


################################################################################
#
# Check for GVars
#
gvars=yes
AC_CHECK_HEADERS(gvars3/instances.h, [], [gvars=no])

AC_CHECK_LIB(GVars3_headless, main)

gvars=no
test "$ac_cv_lib_GVars3_headless_main" == yes && gvars=yes

if test $gvars == yes
then
	gvars=no
	AC_MSG_CHECKING([gvars])
	AC_COMPILE_IFELSE([AC_LANG_SOURCE([
		#include <gvars3/instances.h>
		int main(){GVars3::GUI.LoadFile("");}
	])], [gvars=yes], [gvars=no])

	AC_MSG_RESULT($gvars)
fi

if test $gvars = no
then
	AC_MSG_ERROR([A working copy of GVars3 is required for this program.
Get GVars3 here: http://www.edwardrosten.com/cvd/gvars3.html
])
fi


AC_MSG_CHECKING([TooN support in GVars])
AC_COMPILE_IFELSE([AC_LANG_SOURCE([
	#include <gvars3/instances.h>
	#include <TooN/TooN.h>
	int main(){GVars3::GV3::get<TooN::Vector<> >("", "");}
])], [gvarstoon=yes], [gvarstoon=no])

AC_MSG_RESULT($gvarstoon)
if test $gvarstoon = no
then
	AC_MSG_ERROR([GVars3 has not been compiled with TooN support.
Since TooN is installed, this means you need to recompile GVars3.
])
fi


################################################################################
#
# Check for libcvd
#
cvd=yes
AC_CHECK_HEADERS(cvd/image_io.h, [], [cvd=no])


AC_CHECK_LIB(cvd, main)
cvd=no
cvdlib=no
test "$ac_cv_lib_cvd_main" == yes && cvdlib=yes

if test $cvdlib == yes
then
	AC_MSG_CHECKING(libCVD)
	cvd=no
	AC_LINK_IFELSE([AC_LANG_SOURCE([
		#include <cvd/fast_corner.h>
		#include <cvd/image_io.h>
		int main(){ CVD::Image<CVD::byte> i = CVD::img_load(""); return 0;}
	])], [cvd=yes], [cvd=no])
	AC_MSG_RESULT([$cvd])
fi


if test $cvdlib == yes && test $cvd == no
then
	AC_MSG_CHECKING(If libcvd needs -msse2)
	APPEND(CXXFLAGS,[ -msse2])

	AC_LINK_IFELSE([AC_LANG_SOURCE([
		#include <cvd/fast_corner.h>
		#include <cvd/image_io.h>
		int main(){ CVD::Image<CVD::byte> i = CVD::img_load(""); return 0;}
	])], [cvd=yes], [cvd=no])
	AC_MSG_RESULT($cvd)
fi

if test $cvd = no
then
	AC_MSG_ERROR([A working copy of libCVD is required for this program.
Get libCVD here: http://www.edwardrosten.com/cvd/
])
fi


CHECK_FOR_OPTION(lapack, [
        AC_MSG_CHECKING(if Accelerate framework is needed for LAPACK)
        AC_MSG_RESULT()
        APPEND(LIBS, -framework Accelerate)
        AC_CHECK_FUNC(dgesvd_, [], [a=1])
])


CHECK_FOR_OPTION(gfortran, [
        AC_CHECK_LIB(gfortran, printf, [], [a=1])
])
CHECK_FOR_OPTION(blas, [
        AC_CHECK_LIB(blas, printf, [], [a=1])
])



CHECK_FOR_OPTION(lapack, [
        AC_CHECK_LIB(lapack, dgesvd_, [], [a=1])
])


CHECK_FOR_OPTION(lapack, [
		AC_CHECK_LIB(f2c, main, [], [a=1])
        AC_CHECK_LIB(lapack, dsyev_, [], [a=1])
])


if test "$have_lapack" != yes
then
	AC_MSG_ERROR([Lapack is missing.])
fi


AC_ARG_WITH(imagej, AC_HELP_STRING(--with-imagej=/path/to/imagej/ij.jar, [specify location of ImageJ in order to build plugin]))

if test "$with_imagej" != ""
then
	AC_SUBST(CLASSPATH, [$with_imagej])
fi

java=/usr/lib/jvm/java-6-openjdk/include/

AC_ARG_WITH(jni, AC_HELP_STRING(--with-jni=/path/to/jni, [specify location of JNI headers]))

dnl Test some likely looking places for the JNI header files...

if test "$with_jni" != ""
then
	javas="$with_jni"
else
	jvm="$(readlink -f `which java`  | sed -e's!jre/bin/java!include!')"
	javas="$jvm /usr/lib/jvm/java-7-openjdk/include/ /usr/lib/jvm/java-6-openjdk/include/"
fi

for i in $javas
do
	save_CPPFLAGS="$CPPFLAGS"
	AC_MSG_CHECKING([for jni.h in $i])
	AC_MSG_RESULT([])

	CPPFLAGS="$CPPFLAGS -I$i"

	AC_CHECK_HEADERS([jni.h], [a=0], [a=1])
	unset ac_cv_header_jni_h

	CPPFLAGS="$save_CPPFLAGS"

	if test x$a == x0
	then
		java="$i"
		break
	fi
done

AC_SUBST(JAVAPATH, [$java])

AC_OUTPUT(Makefile)
