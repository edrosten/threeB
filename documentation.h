///@file documentation.h Doxygen documentation bits

/// @defgroup gDebug Useful debugging functions.

/// @defgroup gUtility General utility functions.

/// @defgroup gHMM Generic hidden Markov model solver.
///
/// The notation follows `A tutorial on hidden Markov models and selected applications in speech recognition', Rabiner, 1989
///

/// @defgroup gPlugin Classes related to the ImageJ Plugin

/// @defgroup gStorm Storm classes

// @defgroup gMultiSpotDrift Storm classes specific to multispot processing with drift

/// @defgroup gStormImages Storm imagery classes (basic image processing)

// @defgroup gGraphics Graphical display

// @defgroup gUserInterface General user interface classes

/// @file multispot5.cc Fit spots to the data

/// @file multispot5_headless.cc FitSpots driver for entierly headless (batch) operation

/// @file multispot5_gui.cc FitSpots driver for interactive (GUI) operation and debugging

///@file mt19937.h Mersenne twister interface code

///@file randomc.h 

///@file mersenne.cpp Agner Fogg's Mersenne Twister implementation.

/// @file debug.h  Debugging bits

/// @file debug.cc Debugging bits

/// @file utility.h Utility bits.

/// @file utility.cc Utility bits.

/// @file storm_imagery.cc Code dealing with storm imagery (low level).

/// @file storm_imagery.h Code dealing with storm imagery (low level).

/// @file storm.h Code dealing with storm imagery (high level).

/// @file forward_algorithm.h Contains an implementation fo the forward algorithm.


/// @defgroup gMultiSpot Storm classes specific to multispot processing
// 
// Having completely smooth drift is too memory intensive.
// 
// The drift model is linear, approximated as piecewise constant.  Generally there
// are <i>F</i> frames and these are divided up into <i>S</i> steps. The equation
// for finding the step \e s from the frame \e i is:
// \f[
// i	s = \lfloor \frac{f S}{F} \rfloor.
// \f]
// The inverse is defined to be:
// \f[
// 	f \approx \frac{s + \frac{1}{2}}{S}F
// \f]
//


/** \mainpage 3B Microscopy Analysis

\section sIntro Introduction

This project contains the reference implementation of the 3B microscopy analysis method,
and an ImageJ plugin.
Please refer to <a href="http://www.coxphysics.com/3b">the project website</a> for more information
on the method.

To get started with analysing data, the <a href="http://rsbweb.nih.gov/ij/">ImageJ</a> plugin is
the most suitable piece of software. This can be obtained from the  <a href="http://www.coxphysics.com/3b">the project website</a>.

For more advanced analysis (such as running on a cluster), the 
the program \c multispot5_headless should be used.
This program needs to be run from the commandline.

This project contains the source code for the commandline program and the ImageJ plugin.

\section sStartHere Getting started

\subsection sdata Dataset types and experimental parameters

Bayesian analysis of blinking and bleaching allows data to be extracted from
datasets in which multiple fluorophores are overlapping in each frame.  It can
of course also be used to analyse standard localisation (PALM/STORM) data, but
it must be borne in mind that the low density and high frame number of such
datasets can lead to long runtimes.  Here we briefly discuss the different types
of dataset, related to different applications, that you may wish to analyse
using 3B.

\subsubsection slow Low density PALM/STORM datasets 

These datasets have few fluorophores overlapping
in each image and are at least 10,000 frames long.  As discussed above, they can
be analysed with 3B but their run time makes this a large time investment.  If
you wish to use this approach for performance verification, we would suggest
selecting a small spatial area (around 1.5-3\f$\mu\f$m square).  The algorithm can
also be parallelised by running different sets of a few hundred frames
for the same area
separately.  Even with parallelisation, the large number of
frames will make it time consuming to run.  If you simply want to know what the
structure is like, we would suggest using a method such as QuickPALM, which are
very fast in comparison.

\subsubsection shigh High density fixed cell datasets 

These datasets are of fixed cells but have
multiple fluorophores overlapping in each frame.  You may acquire this type of
dataset if the system you use has a fluorophore which cannot be photoswitched,
or the blinking properties of which cannot be tuned over a wide enough range
using the embedding medium, or if your light source is not powerful enough to
drive most of the fluorophores into the non-emitting state.  In fixed samples
labelled with fluorescent proteins there will in our experience be almost no
blinking present, and 3B will therefore pick up bleaching.  While this can
produce satisfactory results, the localisation is on single event on a high
background and therefore the performance will be significantly degraded.

Some users choose to acquire this type of dense data if they have severe
problems with drift, particularly in the z-direction, as it drastically cuts the
drift over the acquisition time.  However, in the long term it is worth trying
to improve the stability of such systems, since z-drift will impact the accuracy
of all types of localisation measurement.  If you are unsure how badly your
system is affected by z-drift, it is useful to carry out a calibration of the
drift of the microscope using a bead sample before starting acquisition of data
for superresolution.  

\subsubsection slive Live cell datasets 

These datasets are of live cells, generally labelled with
standard fluorescent proteins such as mCherry.  The mounting medium should be
phenol red free, to avoid unnecessary background.  The intensity of the light
source should be selected such that it is high enough to produce blinking but
not strong enough to completely bleach the sample over the time of the
acquisition.  For example, for a standard Xe arc lamp illumination with the full
power of the lamp is generally suitable, but if you are using a powerful laser
source it is recommended to take multiple datasets with different power levels to
determine which is suitable.

The time which is needed to acquire the data necessary for a single
superresolution frame is dependent on the illumination intensity, the speed of
the camera, and the properties of the flurophore.  Many older EMCCD cameras have
a maximum acquisition speed of around 50 frames per second for a small region of
interest, which then gives a limit of 4&nbsp;s to acquire 200 frames.  High
specification EMCCD cameras and SCMOS cameras have much higher acquisition
speeds of up to thousands of frames per second for restricted areas.  However,
in order to maintain the number of photons from each flurophore per frame the
illumination intensity would have to be increased, which is likely to bleach the
sample rapidly, and change the blinking properties of the fluorophore to some
extent.

The acquisition of live cell datasets allows dynamics to be observed.  It should
be noted, however, that if the structures move over the timescale of the
acquisition then that movement will cause blurring in the reconstructed image.


\subsubsection slive Selection of appropriate cell structures for observation 

The 3B algorithm must
be able to pick up the changes in intensity which occur when a fluorophore
switches between an emitting and a non-emitting state.  The accuracy with which
localisation can occur depends, as for other localisation techniques, on the
number of photons from the fluorophore and the background level.  Since 3B is
generally used with a widefield setup, if the sample has a lot of fluorescent
structure out of the plane of focus, the background will be higher and it will
be more difficult to localise.  For thick samples, the background can be reduced
by using TIRF or high angle illumination.



\subsection sanalysis Iterations and run time


In general, determining the number of iterations required for MCMC algorithms
is an unsolved problem. A good general rule is that once the reconstruction
stops changing significantly with increasing iterations, then it is likely that
the reconstruction has converged to a reasonable point.

If you are unsure, then rerun exactly the same area with exactly the same
parameters, but with a different random seed. Note that the ImageJ plugin will
automatically select a different seed each time in the standard interface, but
with the advanced interface or commandline program, the seed must be specified
in the configuration file (see \ref sconfig). If the two results appear
essentially the same, then it is very likely that a sufficient number of
iterations has been reached.

A good general rule is that 200 iterations is sufficient for convergence under
almost all circumstances. For the example usage given in this manual, the
required number of iterations for good convergence requires about 6 hours on a
standard PC (Core i7 at 3GHz).

Convergence may be achieved before 200 iterations, but
terminating the run before convergence can lead to artefacts.  As with any
microscopy method, the experiment and analysis should always be carried out
appropriately to minimise the risk of artefacts.

There are a number of issues which can lead to artefacts:
- Early termination of the algorithm
  - The algorithm builds up a reconstrutcion using a number of random samples.
    If the algorithm is terminated too early, then the reconstruction will be
	dominated by the randomness, and there will not be enough samples for random
	fluctuations to average out.
  - MCMC algorithms may exhibit a property called <a
	href="http://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm">burn
	in</a>. As the 3B algorithm runs, it maintains an estimate of the number of
	spots present in the image. Since it is an MCMC algorithm, this estimate
	will fluctuate around some mean value. However, 3B can add or
	remove spots at a rate of at most 5 per iteration (usually much slower). If
	the algorithm is started with a very bad estimage of the number of spots
	then it may take many iterations for it to approach the mean. During this
	time, the samples drawn will not be representative.
    \n
	3B must therefore run for enough iterations that the later, representative
	iterations will dominate and swamp the earlier ones. If 3B is started using
	a very bad estimate, then this can require substantially more than 200 
	iterations.

- Bright regions close to the boundary
  - The 3B algorithm will not examine any pixels outside the boundary of the
	marked region. If there is a bright region of the image on or near the
	boundary, then 3B will naturally attempt to places fluorophores there. However
	since the point spread function of the microscope is typically several pixels across, the images of fluorophores near
	the boundary will extend across it and so will be missing information, which
	could lead to artefacts.
    \n
	If possible, placing the boundary near a bright region should be avoided. If
	this is not possible, then the reconstruction close to the boundary should
	be ignored. Anything further than the PSF diameter from the boundary is
	unlikely to be affected. Typically this is around 3 pixels.

- Insufficient background regions
  - The 3B algorithm needs to be able to estimate the babkground noise level in
    order to model the image correctly. If there is an insufficient amount of
	background (for instances if the images are too small), then the estimation
	of the noise level may be poor which could lead to artefacts. 
    \n
	The sample data provided is an example of data for which the background can
	be estimated effectively. 

- Image drift or motion
  - The 3B algorithm does not model image motion.
    \n
	If there is significant motion for example due to drift or a live cell
	moving then artefacts may result. The artefacts take the form of streaking
	or smearing in the direction of drift, or structure bunching randomly at one
	end of the drift or the other.
    \n
	Ideally the experiment should be run to minimize drift, but if this is not
	possible, then drift correction software (for example based on tracking
	beads) should be used to correct the drift. For live cell analysis, a tradeoff can be made between spatial and temporal resolution.  If fewer frames are analysed, the temporal resolution is higher but the spatial resolution will be degraded.  Varying the number of frames analysed can also be used to investigate the impact of sample movement in live cells, if this is a concern.

- Incorrect parameters
  - If the parameters (especially the FWHM of the point spread function of the microscope) are set incorrectly then
	3B will not be able to model the image correctly, which may lead to poorer
	resolution or artefacts. The FWHM can be readily determined by taking a
	diffraction limited image of beads. 

- Very high background levels
  - See \ref slive .


\section sUsingPlugin Using the ImageJ Plugin


A tutorial for using the ImageJ plugin is provided under
<code>Plugins>3B>Help</code>

The plugin can operate in two modes, standard and advanced. The standard mode
allows the user to set the microscope PSF and spot size, the starting number of
spots for the analysis and the range of frames.

The plugin also offers an advanced mode of operation which allows much greater
control over 3B. See \ref sconfig for further details.

\section sCMD Using the commandline program

We have provided a set of test data on the website. Download and unpack the zip
file. It will create a new directory called test data with the collowing
contents:
@code 
test_data/AVG_test_data.bmp
test_data/markup1.bmp
img_000000000.fits
img_000000001.fits
img_000000002.fits
...
img_000000299.fits
@endcode

Then run the following command:

@code
./multispot5_headless --save_spots test_data/results.txt --log_ratios test_data/markup1.bmp test_data/img_000000*
@endcode

The program will save the results in the file \c test_data/results.txt . The
program will run indefinitely in the default setup, but you may view the
results at any stage.  There is no well defined stopping point for this type of
algorithm, so it is advisable continuously monitor the resultant image, and
stop the algorithm when the output image is no longer changing with time. After
30 minutes on a fast PC (e.g. Core i7 975), the ring structure which is not
resolved in the widefield image should be clearly visible.
After about 75 mins, the finer details of the structure begin to approach those seen in Fig 2e
in the associated paper.

The ImageJ plugin can load a results file and perform a reconstruction.


Alternatively, you can process the results file further in order to view the results.
Run the following command:
@code
awk '/PASS/{for(i=2; i <=NF; i+=4)print $(i+2), $(i+3)}' test_data/results.txt > test_data/coordinates.txt
@endcode

The file <tt>test_data/coordinates.txt</tt> contains a long list of \f$(x, y)\f$
coordinates, representing possible spots positions. In order to view the
results, load the data into a graph plotting program and create a scatter plot.
NOTE: the axes are in pixel coordinates, so you will have to multiple any
distances by the number of nm/pixel in order to get distances in nm.


\subsection sExplaneExample Example usage explained


\subsubsection ssTestdata Test Data

The 300 TIFF files in the test directory correspond to the data used for Fig. 2
in the paper. Please refer to the paper for details on how the data was
obtained.

The file <tt>AVG_test_data.bmp</tt> is a Z projection made using 
<a href="http://rsbweb.nih.gov/ij/">ImageJ</a>.

The file <tt>markup1.bmp</tt> is a mask indicating which area of the image to
analyse. All perfectly black pixels are ignored, ecerything else is analysed. If
you overlay <tt>markup1.bmp</tt> and <tt>AVG_test_data.bmp</tt> you can see
which area the markup corresponds to. The markup file was created using 
<a href="http://www.gimp.org">the GIMP</a>.

\subsubsection ssRunning Running the program

The general form for running the program is:

@code
./multispot5_headless [ --variable1 value1 [ --variable2 value 2 [ ... ] ] ] image1 image2 ... 
@endcode

so the example sets ths variable \c save_spots to \c test_data/results.txt and
the variable \c log_ratios to \c test_data/markup1.bmp.  The remaining
arguments is the list of files to be analysed.

The program gets the markup in the filename given in the \c log_ratios variable
(yes, the choice of name is very strange, and corresponds to a very old phase of
development). The more sanely named variable \c save_spots is the filename in
which the output is to be saved.

The program actually has a large number of variables which must be set. Most of
them you probably don't want to change, but some of them you will want to
change. The default values for these variables are stored in \c multispot5.cfg
The format of this file should be mostly self explanatory. Everything after a \c
// is a comment and is ignored. See \ref sconfig for further details.

You will probably want to change:
<ul>
<li> \c blur.mu

This is the prior over spot size, which is how the pixel size and microscope
FWHM are represented. Some example values for a FWHM of 300nm/pix at 160 and
100 nm/pix and for a FWHM of 270nm at 79nm per pixel. 

If you have significantly largre or smaller pixels, the performance may be
degraded.

<li> \c placement.uniform.num_spots

This is the initial number of spots to be placed down. Eventually, the algorithm
will converge to a reasonable number of spots, even if this value is far off.
The default value (15) is appropriate given the small area and dimness of the
sample data. You will want to increase this number for larger areas of markup
and relatively brighter regions.

If this number is more than 1000, then the algorithm will run very slowly and
may take several days.

</ul>

Note that variables specified on the commandline override all variables in the
configuration file.

The program can read FITS, BMP, PPM and PGM images. Depending on how it
has been compiled, it can also read TIFF, PNG and JPEG images.
The program cannot work on multi-image TIFF files. ImageJ can be used to split a
multi-image TIFF into a collection of single image files. All the images loaded
must be the same size.

\subsubsection ssExtract Extracting and visualising the data

The output file (in this case  \c test_data/results.txt ) containing the results
is in a format unsuitable for plotting directly, and must be extracted. The
reason for this is that the output file contains enough information to
seamlessly continue long runs which have been interrupted. The provided AWK
program extracts the coordinates of the spots over all iterations and puts them
in \c coordinates.txt.

Alternatively, the data can ve visualised using the plugin under the menu
<code>Plugins>3B>Open 3B run</code>

\section sconfig The configuration file and advanced settings

The 3B system has a large number of parameters which control its behaviour.
These are controled via the configuration file for the commandlie program or
via the ``Advanced'' option for the plugin. The ``Advanced'' option essentially
allows you to supply a fully custom configuration file.

The sample configuration file is given below along with explanations of all
parameters. In the program, most of these parameters are used by the 
FitSpots class.

\include jar/multispot5.cfg



\page sComp Compiling the programs

In order to compile the project, you will need to download and install the
following libraries:
<ul>
<li> TooN http://www.edwardrosten.com/cvd/toon.html
<li> libcvd http://www.edwardrosten.com/cvd/
<li> gvars3 http://www.edwardrosten.com/cvd/gvars3.html
</ul>
The program is portable and is well tested under Linux and OSX. It will also
compile under Windows using cygwin or MinGW.

The program can be built using the usual method for compiling under Linux:
@code
./configure && make 
@endcode

\section sPlugin Compiling the ImageJ plugin

The plugin is provided pre-compiled from the project website.

There are two ways of building the plugin, manual and automatic. If you want to
make changes to the plugin, then use manual building. If you want to
automatically build the plugin for several platforms, then use automatic
building.

\subsection sManual Manual

The basic build instructions are the same as for the commandline program.
You will also need the JDK (Java Development Kit) and  ImageJ installed.

You will have to locate where your system has installed the JDK. If it is 
not in /usr/lib/jvm/java-6-openjdk/include, you will have to specify the path:

First configure the system:
@code
./configure --with-imagej=/path/to/ImageJ/ij.jar --with-jni=/path/to/jdk/include
@endcode

You will also need to make sure that the JDK programs (javac, havah, etc) are in
your path. The configure script will attempt to detect the location of the JNI headers.
If it fails, you will need to specify \c --with-jni=/path/to/jdk/include


To build the JAVA part:
@code
make three_B.jar
@endcode


To build the plugin (on Linux):
@code
make libthreeB_jni.so DYNAMIC_PLUGIN=1
@endcode
Note that if you do not specify \c DYNAMIC_PLUGIN, then the makefile will try to
build a plugin with some dependencies statically linked in which will almost
certainly fail unless you have set the system up to support  such an operation.

On MinGW:
@code
make threeB_jni.dll
@endcode

Now copy three_B.jar and libthreeB_jni.so into your ImageJ plugins directory.


\subsection sAutoBuild Automatic Build

The automatic build method is very slow and is designed to be able to repeatably 
build plugins for 32 and 64 bit Linux and Windows. It is also designed to build
the plugin with as many static dependencies as possible so that only a single
DLL/so needs to be shipped per system.

The script operates by building a temporary install of Ubuntu 10.04 LTS, and
using that to compile all variants of the plugin. 

You will need a Debian based system (or a system on which the command \c
debootstrap works) and root access.

Tha automatic build system makes use of cLAPACK, rather than LAPACK as the
LAPACK part is not speed critical and it is easier to build CLAPACK without
additional external dependencies.

To build, run the following commands:
@code
#First make a tar.gz of the source code
bash make_dist.sh

#Now execute the automatic build process
bash build_plugin.sh
@endcode

The build takes a long time, and you should probably edit \c build_plugin.sh to
point the installer at an Ubuntu mirror somewhere near to where you are.

At the end of the build, the script will print out a directory name like:
@code
dist-123908
@endcode

A fresh copy of the plugin DLL and shared object will be present in that
directory named.


*/
