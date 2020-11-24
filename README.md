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


Building
--------

If you wish to run 3B on a cluster, then you need to either build a static
executable or take care to copy all the dependencies. A script has been
provided To build for all platforms. You will need a Linux machine. On ubuntu
or debian, run:

```
sudo apt install debootstrap
sudo bash build_plugin.sh
```

Note, you need to have the `debootstrap` program available, so you will need
this from another source if you are not on Ubuntu.


sudo access is needed since the script runs debootstrap which is used to make
some clean, consistent debian environments for a repeatable build. This is
instead of using a system like Docker which was not widespread and mature when
this was originally written.


Building
--------

The usual ./configure && make can be used to build the software.
The software depends on:

TooN http://www.edwardrosten.com/cvd/toon.html
libcvd http://www.edwardrosten.com/cvd/
gvars3 http://www.edwardrosten.com/cvd/gvars3.html

You currently need old versions (newer ones will work, but not too new)

http://www.edwardrosten.com/cvd/libcvd-20121025.tar.gz
http://www.edwardrosten.com/cvd/gvars-3.0.tar.gz
http://www.edwardrosten.com/cvd/TooN-2.0.tar.gz


License
-------

The software may be redistributed under the terms of the GPL 3.0

