//The functiins from storm_imagery.h which don't make use of img_load
//Separated to allow JNI plugin to build more easily.
#include "debug.h"
#include "utility.h"
#include "storm_imagery.h"
#include <gvars3/instances.h>
#include <cvd/convolution.h>

using namespace CVD;
using namespace GVars3;
using namespace std;

/**Do the initial preprocessing on a loaded image. Currently
this is a high pass filter to make the resultimg images zero mean.

The filter is controlled with the \c preprocess.lpf and \c preprocess.skip Gvars

See also load_and_preprocess_images()

@param im  images
@return preprocessed images.
@ingroup gStormImages
*/
Image<float> preprocess_image(const Image<float>& im)
{
	float wide = GV3::get<float>("preprocess.lpf", 0., -1);
	bool p = GV3::get<bool>("preprocess.skip", 0, -1);

	//Highpass filter the images using blur and subtract
	if(!p)
	{
		Image<float> f(im.size(), 0), fwide(im.size(), 0);
		convolveGaussian_fir(im, fwide, wide);
		for(int r=1; r < im.size().y-1; r++)
			for(int c=1; c < im.size().x-1; c++)
				f[r][c] = im[r][c] - fwide[r][c];

		return f;
	}
	else
		return im.copy_from_me();
}

/**Find the mean and variance of a stack of images
@param images Image stack
@return (mean, variance)
@ingroup gStormImages
*/
pair<float, float> mean_and_variance(const vector<Image<float> >& images)
{
	assert_same_size(images);

	double sum=0, sum2=0, area=0;

	for(unsigned int i=0; i < images.size(); i++)
	{
		area += images[i].size().area();
		for(int r=0; r < images[i].size().y; r++)
			for(int c=0; c < images[i].size().x; c++)
			{
				sum += images[i][r][c];
				sum2 += sq(images[i][r][c]);
			}
	}

	sum /= area;
	sum2 /= area;
	return make_pair(sum, sum2 - sq(sum));
}


