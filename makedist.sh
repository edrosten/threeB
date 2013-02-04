echo "Makedist! $1"

autoconf
./configure
make clean
ccfiles=`make -nB all | awk '{for(i=1;i<=NF; i++)if($(i)~/\.(cc)|(cpp)$/)print $(i)}'`
hfiles=`make -nB all | awk '{for(i=1;i<=NF; i++)if($(i)~/\.(cc)|(cpp)$/)print $(i)}' | xargs makedepend -f - 2> /dev/null  | awk '{for(i=1;i<=NF; i++)if($(i)~/\.hh?$/)print $(i)}' | sort -u | grep -v usr`

dir=$1

if ! mkdir $dir
then
	echo "Error $dir exists."
	exit 1
fi


for i in $ccfiles $hfiles 
do
	cat header.txt $i > $dir/$i
done

cp *.java $dir
cp ThreeBRunner.h $dir
cp -r pstreams $dir
cp three_B.jar $dir
cp -r jar $dir
cp -r tag $dir
cp LICENSE README $dir
cp documentation.h Doxyfile storm.sty $dir
cp multispot5.cfg spot_viewer.cfg $dir
cp configure.ac configure Makefile.in $dir

cd $dir
#doxygen
cd ..
tar -cvzf $dir.tar.gz $dir/
tar -cvjf $dir.tar.bz2 $dir/
tar -cvJf $dir.tar.xz $dir/
zip -qr $dir.zip $dir


