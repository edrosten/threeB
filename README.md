ThreeB
------

This is the software for 3B microscopy analysis.

For information on using the software, please refer to the documentation (in the
html/ directory) or http://www.coxphysics.com/3b/index.html


If you want to try or test the software, then the recomended method is to use a
release of the ImageJ plugin (which supports Windows and Linux). If you want
more control or to run on a cluster, then the recommended method is to use a
precompiled binary from the most recent release.

If you download a precompiled executable for a cluster, you will still need the
`multispot5.cfg` configuration file from this repository. 

Relseaes
--------

http://www.coxphysics.com/3b/3B_releases/


Building the static executables
-------------------------------

The easiest and fastest way is using docker, so you need to get docker. Then go into
your source directory and run:

```
docker run -v $PWD:/home/build edrosten/threeb-static-build-env ./configure
docker run -v $PWD:/home/build edrosten/threeb-static-build-env make multispot5_static -j 8
```


Building the plugin (and static executables)
--------------------------------------------

This is the "classic" way of building 3B, and I haven't ported this over to docker
yet. It does more or less the same thing, but predates docker by about 6 years or 
so, so you need sudo access.

You will need a Linux machine. On ubuntu or debian, run:

```
sudo apt install debootstrap
sudo bash build_plugin.sh
```

sudo access is needed since the script runs debootstrap which is used to make
some clean, consistent debian environments for a repeatable build.


Building completely by hand
---------------------------

Don't do this unless you are really very confident and 100% sure you need it.

Building the is done the usual way (`./configure && make`), but that's a little
tricky these days because (1) building static binaries is quite hard and (2) I 
wrote this in 2009 and released it in 2011, so it doesn't currently build on 
new systems with up to date library versions.

The usual ./configure && make can be used to build the software.
The software depends on:

TooN http://www.edwardrosten.com/cvd/toon.html
libcvd http://www.edwardrosten.com/cvd/
gvars3 http://www.edwardrosten.com/cvd/gvars3.html

You currently need old versions (newer ones will work, but not too new)

http://www.edwardrosten.com/cvd/libcvd-20121025.tar.gz
http://www.edwardrosten.com/cvd/gvars-3.0.tar.gz
http://www.edwardrosten.com/cvd/TooN-2.0.tar.gz

You'll likely need an old compiler as well. The static executables will build with
gcc-4.4.3. If you want to build the plugin, you'll need a newer MinGW environment.
Building for MinGW is hard and requires patching then building dependencies.

License
-------

The software may be redistributed under the terms of the GPL 3.0

