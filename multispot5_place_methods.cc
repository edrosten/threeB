#include <tag/printf.h>
#undef make_tuple

#include <climits>
#include <algorithm>
#include <tr1/tuple>
#include <cvd/image.h>
#include <cvd/vector_image_ref.h>
#include <cvd/byte.h>
#include <gvars3/instances.h>
#include <TooN/TooN.h>

#include "multispot5.h"
#include "forward_algorithm.h"
#include "storm.h"
#include "debug.h"

using namespace std;
using namespace CVD;
using namespace GVars3;
using namespace TooN;
using namespace std::tr1;


tuple<StateParameters, double, double> static get_defaults(const vector<ImageRef>& pixels)
{
	const double variance = 1; // it should be

	//To scale the X axis of a log-normal distribution, only
	//the mu parameter needs to be changed...
	const double intensity_mu = GV3::get<double>("intensity.rel_mu", 0., -1)  + log(sqrt(variance));
	const double intensity_sigma = GV3::get<double>("intensity.rel_sigma", 0., -1);
	const double blur_mu = GV3::get<double>("blur.mu", 0., -1);
	const double blur_sigma = GV3::get<double>("blur.sigma", 0., -1);


	const double intensity  = log_normal_mode(intensity_mu, intensity_sigma);
	const double blur =  log_normal_mode(blur_mu, blur_sigma);

	StateParameters p;

	//Initialize the MT19937 RNG from a seed.
	p.rng = shared_ptr<MT19937>(new MT19937);
	p.rng->simple_seed(GV3::get<int>("seed", 0, 1));

	//Start at the beginning
	p.pass=0;
	p.iteration=0;

	//Remember which pixels are in use
	p.pixels = pixels;

	return make_tuple(p, intensity, blur);
}

//Function for accumulating an integer (counting)
static void acc(int & i, const Vector<2>&)
{
	i++;
}

//Function for accumulating a container (inserting)
static void acc(vector<Vector<2> >& i, const Vector<2> &r)
{
	i.push_back(r);
}

//Function for placing spots over a hexagonal grid 
template<class Ret>
Ret place_spots(double sp, Vector<2> centre, double radius, const Image<bool>& mask)
{
	double angle=M_PI/180 * 6;
	Vector<2> a_axis = SO2<>(angle) * makeVector(1,0);
	Vector<2> b_axis = SO2<>(M_PI/3) * a_axis;
	
	Ret num = Ret();
	//The range is +- 2*r / sqrt(3) on each axis.
	//To prove:
	//Draw a curcle with a horizontal line through it,
	//brushing the top and brushing the bottom (the a axis).
	//
	//Now draw the same three lines, but rotated round by 60 degrees.
	//(the b axis).
	//
	//A bunch of 30/60/90 triangles with an opposite length of r
	//are formed. We need the hypotenuse which is 2r/sqrt(3).

	int n = (int) ceil(2*radius/sqrt(3) / sp);
	for(int a=-n; a <= n; a++)
		for(int b=-n; b <= n; b++)
		{
			Vector<2> pos = centre + a*sp*a_axis + b*sp*b_axis;
			ImageRef p = ir(pos+makeVector(.5, .5));

			if(mask.in_image(p) && mask[p])
				acc(num, pos);
		}
	
	return num;
}

vector<Vector<2> > find_spacing(int target, const Image<bool>& mask)
{	
	//First a bounding circle is required. The circle need not be tight.

	//Compute the approximate bounding circle
	//First compute the centroid
	Vector<2> centre = Zeros;
	int count=0;
	for(int y=0; y < mask.size().y; y++)
		for(int x=0; x < mask.size().x; x++)
			if(mask[y][x])
			{
				centre += makeVector(x, y);
				count ++;
			}
	centre /= count;
	
	double r2 = 0;
	//Now compute the radius
	for(int y=0; y < mask.size().y; y++)
		for(int x=0; x < mask.size().x; x++)
			if(mask[y][x])
				r2 = max(r2, norm_sq(makeVector(x,y) - centre));
	double radius = r2;

	//Perform a binary search to find an appropriate spacing. The function
	//is not monotonic, so the spacing is not guaranteed to be unique.
	double close = 0;
	int large_num = INT_MAX;

	double far = sqrt(mask.size().mag_squared());
	int small_num = place_spots<int>(far, centre, radius, mask);
	
	if(target > small_num)
		while(small_num != large_num && far - close > 1e-6)
		{
			double mid = (close + far)/2;
			int mid_num = place_spots<int>(mid, centre, radius, mask);
			
			if(mid_num > target)
			{
				large_num = mid_num;
				close = mid;
			}
			else
			{
				small_num = mid_num;
				far = mid;
			}
		}
	
	//Pick the best, in case the algorithm terminated due to 
	//a too small disparity.
	double spacing;
	if(large_num - target < target - small_num)
		spacing = close;
	else
		spacing = far;
	
	//Use small_num or close as the spacing
	return place_spots<vector<Vector<2> > >(spacing, centre, radius, mask);
}


StateParameters place_spots_uniform(int num_spots, const vector<ImageRef>& pixels, const ImageRef& size)
{
	Image<bool> pix(size);
	pix.fill(false);
	for(unsigned int i=0; i < pixels.size(); i++)
		pix[pixels[i]] = true;


	vector<Vector<2> > spot_pos = find_spacing(num_spots, pix);


	StateParameters p;
	double intensity_mode, blur_mode;
	tie(p, intensity_mode, blur_mode) = get_defaults(pixels);

	//Create all the spots
	for(unsigned int i=0; i < spot_pos.size(); i++)
		p.spots.push_back(makeVector(intensity_mode, blur_mode, spot_pos[i][0], spot_pos[i][1]));

	
	return p;
}

StateParameters place_spots_intensity_sampled(int num_spots, const vector<ImageRef>& pixels, const vector<Image<float> >& ims)
{
	assert_same_size(ims);

	StateParameters p;
	double intensity_mode, blur_mode;
	tie(p, intensity_mode, blur_mode) = get_defaults(pixels);
	
	//Assume that the average intensity is a unnormalized probability distribution
	//and use rejection sampling to sample uniformly over it!
	
	//No need to scale...
	vector<float> intensities(pixels.size(),0);
	for(unsigned int i=0; i < pixels.size(); i++)
		for(unsigned int j=0; j < ims.size(); j++)
			intensities[i] += ims[j][pixels[i]];

	double max_intensity = *max_element(intensities.begin(), intensities.end());

	if(max_intensity < 0)
		return p;

	MT19937& rng = *(p.rng);

	while((int)p.spots.size() < num_spots)
	{
		int element = floor(rng() * pixels.size());
		double y = rng() * max_intensity;

		if(y <= intensities[element])
			p.spots.push_back(makeVector(intensity_mode, blur_mode, pixels[element].x + rng()-.5, pixels[element].y + rng()-.5));
	}

	return p;
}
