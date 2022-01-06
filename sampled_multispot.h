#ifndef SPOT_WITH_BACKGROUND_H
#define SPOT_WITH_BACKGROUND_H

#include <vector>
#include <cvd/image_ref.h>
#include <tuple>
#include <TooN/TooN.h>

#include "drift.h"

typedef char State;

namespace SampledMultispot
{

using namespace std;
using namespace CVD;
using namespace TooN;



//The changes to SpotWithBackground to operate with/without drift and masking are
//irritatingly pervasive. I can't think of any other way of doing it.
#define SWBG_NAME SpotWithBackground
#include "spot_with_background.hh"


#define SWBG_NAME SpotWithBackgroundMasked
#define SWBG_HAVE_MASK
#include "spot_with_background.hh"

//#define SWBG_NAME SpotWithBackgroundDrift
//#define SWBG_HAVE_DRIFT
//#include "spot_with_background.hh"

//#undef SWBG_NAME
//#define SWBG_NAME SpotWithBackgroundDriftMasked
//#define SWBG_HAVE_DRIFT
//#define SWBG_HAVE_MASK
//#include "spot_with_background.hh"




inline double intensity(double i)
{
	return i;
}

inline double intensity(const pair<double, Vector<4> >& i)
{
	return i.first;
}


//Add and remove a spot over the entire region
template<class T>
void remove_spot(vector<vector<double> >& current_sample_intensities, const vector<T>& spot_intensities, const vector<State>& spot_sample)
{
	for(unsigned int frame=0; frame < current_sample_intensities.size(); frame++)
		if(spot_sample[frame] == 0) //Spot is on, so remove it
			for(unsigned int p=0; p < spot_intensities.size(); p++)
				current_sample_intensities[frame][p] -= intensity(spot_intensities[p]);
}

template<class T>
void add_spot(vector<vector<double> >& current_sample_intensities, const vector<T>& spot_intensities, const vector<State>& spot_sample)
{
	for(unsigned int frame=0; frame < current_sample_intensities.size(); frame++)
		if(spot_sample[frame] == 0) //Spot is on, so add it
			for(unsigned int p=0; p < spot_intensities.size(); p++)
				current_sample_intensities[frame][p] += intensity(spot_intensities[p]);
}


//Add and remove a spot only over a mask. 
template<class T>
void remove_spot(vector<vector<double> >& current_sample_intensities, const vector<T>& spot_intensities, const vector<State>& spot_sample, const vector<int>& mask)
{
	for(unsigned int frame=0; frame < current_sample_intensities.size(); frame++)
		if(spot_sample[frame] == 0) //Spot is on, so remove it
			for(unsigned int p=0; p < mask.size(); p++)
				current_sample_intensities[frame][mask[p]] -= intensity(spot_intensities[mask[p]]);
}

template<class T>
void add_spot(vector<vector<double> >& current_sample_intensities, const vector<T>& spot_intensities, const vector<State>& spot_sample, const vector<int>& mask)
{
	for(unsigned int frame=0; frame < current_sample_intensities.size(); frame++)
		if(spot_sample[frame] == 0) //Spot is on, so add it
			for(unsigned int p=0; p < mask.size(); p++)
				current_sample_intensities[frame][mask[p]] += intensity(spot_intensities[mask[p]]);
}


//Add and remove a drifty spot only over a mask.
template<class T>
void remove_spot(vector<vector<double> >& current_sample_intensities, const vector<vector<T> > & spot_intensities, const vector<State>& spot_sample, const vector<int>& mask)
{
	const int steps = spot_intensities.size();
	const int frames = current_sample_intensities.size();

	for(int frame=0; frame < frames; frame++)
	{
		int s = frame * steps / frames;

		if(spot_sample[frame] == 0) //Spot is on, so remove it
			for(unsigned int p=0; p < mask.size(); p++)
				current_sample_intensities[frame][mask[p]] -= intensity(spot_intensities[s][mask[p]]);
	}
}

template<class T>
void add_spot(vector<vector<double> >& current_sample_intensities, const vector<vector<T> >& spot_intensities, const vector<State>& spot_sample, const vector<int>& mask)
{
	const int steps = spot_intensities.size();
	const int frames = current_sample_intensities.size();

	for(int frame=0; frame < frames; frame++)
	{
		int s = frame * steps / frames;

		if(spot_sample[frame] == 0) //Spot is on, so add it
			for(unsigned int p=0; p < mask.size(); p++)
				current_sample_intensities[frame][mask[p]] += intensity(spot_intensities[s][mask[p]]);
	}
}



//Compute the spot intensity for a given spot at each pixel
inline vector<double> compute_spot_intensity(const vector<ImageRef>& pixels, const Vector<4>& params)
{
	vector<double> intensities(pixels.size());

	for(unsigned int i=0; i < pixels.size(); i++)
		intensities[i] = spot_shape(vec(pixels[i]), params);

	return intensities;
}

//Compute the spot intensity derivatives for a given spot at each pixel
inline vector<pair<double, Vector<4> > > compute_spot_intensity_derivatives(const vector<ImageRef>& pixels, const Vector<4>& params)
{
	vector<pair<double, Vector<4> > > derivatives(pixels.size());

	for(unsigned int i=0; i < pixels.size(); i++)
		derivatives[i] = spot_shape_diff_position(vec(pixels[i]), params);
	return derivatives;
}

inline vector<tuple<double, Vector<4>, Matrix<4> > > compute_spot_intensity_hessian(const vector<ImageRef>& pixels, const Vector<4>& params)
{
	vector<tuple<double, Vector<4>, Matrix<4> > > hessian(pixels.size());

	for(unsigned int i=0; i < pixels.size(); i++)
		hessian[i] = spot_shape_hess_position(vec(pixels[i]), params);
	return hessian;
}


/**
Create a sequence of integers. These can be used as observations
in an observation class by forward_algorithm() and etc.
@param n Length of sequence
@ingroup gUtility
*/
inline vector<int> sequence(int n)
{
	vector<int> v;
	for(int i=0; i < n; i++)
		v.push_back(i);
	return v;
}

/*struct RndGrand48
{
	double operator()()
	{
		return drand48();
	}
};*/

///Draw samples from the spot states given the spots positions and some data.
///Variable naming matches that in FitSpots.
///@ingroup gStorm
class GibbsSampler
{
	const vector<vector<double> >& pixel_intensities;
	const vector<vector<double> >& spot_intensities;
	const vector<Vector<4> > spots;
	const Matrix<3> A; 
	const Vector<3> pi;
	const double base_variance;
	double variance;

	const int sample_iterations;
	const int num_frames, num_pixels;
	const vector<int> O;
  
	vector<vector<State> > current_sample;
	vector<vector<double> > current_sample_intensities;

	public:
	
	GibbsSampler(const vector<vector<double> >& pixel_intensities_,
	             const vector<vector<double> >& spot_intensities_,
	             const vector<Vector<4> >& spots_,
	             const Matrix<3> A_,
	             const Vector<3> pi_,
	             double variance_,
	             int sample_iterations_)
	:pixel_intensities(pixel_intensities_), //pixel_intensities: [frame][pixels]
	 spot_intensities(spot_intensities_),   //spot_intensities: [spot][pixel]
	 spots(spots_),
	 A(A_),
	 pi(pi_),
	 base_variance(variance_),
	 variance(variance_),
	 sample_iterations(sample_iterations_),
	 num_frames(pixel_intensities.size()),
	 num_pixels(pixel_intensities[0].size()),
		//Observations vector. As usual for this application, the observations are just
		//numbered integers which refer to data held elsewhere.
	 O(sequence(num_frames)),
		//Start all spots OFF, so the intensity is 0. OFF is 1 or 2, not 0!!!
		//sample_list: [sample][spot][frame]: list of samples drawn using Gibbs sampling
	 current_sample(spots.size(), vector<State>(num_frames, 2)), //current sample [spot][frame]
		//pixel intensities assosciated with the current sample [frame][pixel]
	 current_sample_intensities(num_frames, vector<double>(num_pixels))
	{
		//Check a bunch of stuff
		assert_same_size(pixel_intensities);
		assert_same_size(spot_intensities);

	}

	///Update the noide variance. Used for adding thermal noise.
	///@param v noise variance.
	void set_variance(double v)
	{
		variance = v;
	}


	///Reset the gibbs sampler oro the initial state (all spots off)
	void reset()
	{
		vector<State> off(num_frames, 2);
		fill(current_sample.begin(), current_sample.end(), off);

		vector<double> black(num_pixels);
		fill(current_sample_intensities.begin(), current_sample_intensities.end(), black);
		variance = base_variance;
	}
	
	///Get the next sample
	///@param rng Random number generator
	template<class T> void next(T& rng)
	{
		for(int j=0; j < sample_iterations; j++)
			for(int k=0; k < (int) spots.size(); k++)
			{
				//Subtract off the spot we're interested in.
				remove_spot(current_sample_intensities, spot_intensities[k], current_sample[k]);

				//Now current_sample_intensities is the image value for every spot in every frame,
				//except the current spot, which is always set to off. This allows us to add it in 
				//easily.
				SpotWithBackground B(current_sample_intensities, spot_intensities[k], pixel_intensities, variance);
				vector<array<double, 3> > delta = forward_algorithm_delta(A, pi, B, O);
				current_sample[k] = backward_sampling<3,State, T>(A, delta, rng);

				//Put the newly sampled spot in
				add_spot(current_sample_intensities, spot_intensities[k], current_sample[k]);
			}
	}
	/*void next()
	{
		RngDrand48 rng;
		next(rng);
	}*/	

	///Retrieve the current sample
	const vector<vector<State> >& sample() const
	{
		return current_sample;
	}
	///Retrieve the intensities for the current sample
	const vector<vector<double> >& sample_intensities() const
	{
		return current_sample_intensities;
	}

};

///Gibbs sampling class which masks spots to reduce computation.
///
///This draws samples from, the spot states given the spots positions and some data. It is
///very similar to GibbsSampler, except that it only computes probabilities in a mask around each spot
///to save on computation. Variable naming matches that in FitSpots.
///@ingroup gStorm
class GibbsSampler2
{
	const vector<vector<double> >& pixel_intensities;
	const vector<vector<double> >& spot_intensities;
	const vector<Vector<4> > spots;
	const std::vector<std::vector<int> >& spot_pixels;
	const Matrix<3> A; 
	const Vector<3> pi;
	const double base_variance;
	double variance;

	const int sample_iterations;
	const int num_frames, num_pixels;
	const vector<int> O;
  
	vector<vector<State> > current_sample;
	vector<vector<double> > current_sample_intensities;

	vector<double> cutout_spot_intensities;
	vector<vector<double> > cutout_pixel_intensities;
	vector<vector<double> > cutout_current_sample_intensities;

	public:
	
	GibbsSampler2(const vector<vector<double> >& pixel_intensities_,
	             const vector<vector<double> >& spot_intensities_,
	             const vector<Vector<4> >& spots_,
	             const vector<vector<int> >& spot_pixels_,
	             const Matrix<3> A_,
	             const Vector<3> pi_,
	             double variance_,
	             int sample_iterations_)
	:pixel_intensities(pixel_intensities_), //pixel_intensities: [frame][pixels]
	 spot_intensities(spot_intensities_),   //spot_intensities: [spot][pixel]
	 spots(spots_),
	 spot_pixels(spot_pixels_),
	 A(A_),
	 pi(pi_),
	 base_variance(variance_),
	 variance(variance_),
	 sample_iterations(sample_iterations_),
	 num_frames(pixel_intensities.size()),
	 num_pixels(pixel_intensities[0].size()),
		//Observations vector. As usual for this application, the observations are just
		//numbered integers which refer to data held elsewhere.
	 O(sequence(num_frames)),
		//Start all spots OFF, so the intensity is 0. OFF is 1 or 2, not 0!!!
		//sample_list: [sample][spot][frame]: list of samples drawn using Gibbs sampling
	 current_sample(spots.size(), vector<State>(num_frames, 2)), //current sample [spot][frame]
		//pixel intensities assosciated with the current sample [frame][pixel]
	 current_sample_intensities(num_frames, vector<double>(num_pixels)),
	 cutout_pixel_intensities(num_frames),
	 cutout_current_sample_intensities(num_frames)
	{
		//Check a bunch of stuff
		assert_same_size(pixel_intensities);
		assert_same_size(spot_intensities);
	}

	///Update the noide variance. Used for adding thermal noise.
	///@param v noise variance.
	void set_variance(double v)
	{
		variance = v;
	}


	///Reset the gibbs sampler oro the initial state (all spots off)
	void reset()
	{
		vector<State> off(num_frames, 2);
		fill(current_sample.begin(), current_sample.end(), off);

		vector<double> black(num_pixels);
		fill(current_sample_intensities.begin(), current_sample_intensities.end(), black);
		variance = base_variance;
	}
	
	///Get the next sample
	///@param rng Random number generator
	template<class T> void next(T& rng)
	{
	
//double remove=0;
//double cut=0;
//double swb=0;
//double ff_masked=0;
//double bs=0;
//double add=0;
//cvd_timer t;
	std::vector<array<double, 3> > delta3;
		for(int j=0; j < sample_iterations; j++)
			for(int k=0; k < (int) spots.size(); k++)
			{
//t.reset();
				//Subtract off the spot we're interested in.
				remove_spot(current_sample_intensities, spot_intensities[k], current_sample[k], spot_pixels[k]);
//remove+=t.reset();
/*
				//Cut out
				//spot
				cutout_spot_intensities.resize(spot_pixels[k].size());
				for(unsigned int i=0; i < spot_pixels[k].size(); i++)
					cutout_spot_intensities[i] = spot_intensities[k][spot_pixels[k][i]];

				//others
				for(int f=0; f < num_frames; f++)
				{
					cutout_current_sample_intensities[f].resize(spot_pixels[k].size());
					cutout_pixel_intensities[f].resize(spot_pixels[k].size());
					for(unsigned int i=0; i < spot_pixels[k].size();i++)
					{
						cutout_current_sample_intensities[f][i] = current_sample_intensities[f][spot_pixels[k][i]];
						cutout_pixel_intensities[f][i] = pixel_intensities[f][spot_pixels[k][i]];
					}
				}*/
//cut += t.reset();
				//Now current_sample_intensities is the image value for every spot in every frame,
				//except the current spot, which is always set to off. This allows us to add it in 
				//easily.

//				SpotWithBackground B(current_sample_intensities, spot_intensities[k], pixel_intensities, variance);
//				vector<array<double, 3> > delta = forward_algorithm_delta(A, pi, B, O);

//ff+=t.reset();

//				SpotWithBackground B2(cutout_current_sample_intensities, cutout_spot_intensities, cutout_pixel_intensities, variance);
//				std::vector<array<double, 3> > delta2 = forward_algorithm_delta(A, pi, B2, O);
//ff_cut+=t.reset();

				SpotWithBackgroundMasked B3(current_sample_intensities, spot_intensities[k], pixel_intensities, variance, spot_pixels[k]);
//swb += t.reset();
				forward_algorithm_delta2<3>(A, pi, B3, O, delta3);
//f_masked+=t.reset();
				/*for(unsigned int i=0; i < delta.size(); i++)
				{
					cout.precision(20);
					cout.setf(cout.scientific);
					std::cout << delta[i][0] << " " << delta[i][1] << " " <<delta[i][2] << std::endl;
					std::cout << delta2[i][0] << " " << delta2[i][1] << " " <<delta2[i][2] << std::endl;
					cout << endl;
				}
				std::exit(1);

*/

				current_sample[k] = backward_sampling<3,State, T>(A, delta3, rng);
//bs += t.reset();
				//Put the newly sampled spot in
				add_spot(current_sample_intensities, spot_intensities[k], current_sample[k], spot_pixels[k]);
//add += t.reset();
			}
//		cout << "remove=" <<remove << " cut=" << cut << " swb=" << swb<< " ff_mask=" << ff_masked << " bs=" <<bs << " add="<<add << endl;
	}
/*	void next()
	{
		RngDrand48 rng;
		next(rng);
	}	
*/	
	///Retrieve the current sample
	const vector<vector<State> >& sample() const
	{
		return current_sample;
	}
	///Retrieve the intensities for the current sample
	const vector<vector<double> >& sample_intensities() const
	{
		return current_sample_intensities;
	}

};

#if 0
///Gibbs sampling class
class GibbsSampler3
{
	const vector<vector<double> >& pixel_intensities;
	const vector<vector<vector<double> > >& spot_intensities;
	const vector<Vector<4> > spots;
	const std::vector<std::vector<int> >& spot_pixels;
	const Matrix<3> A; 
	const Vector<3> pi;
	const double base_variance;
	double variance;

	const int sample_iterations;
	const int num_frames, num_pixels;
	const vector<int> O;
  
	vector<vector<State> > current_sample;
	vector<vector<double> > current_sample_intensities;

	public:
	
	GibbsSampler3(const vector<vector<double> >& pixel_intensities_,
	             const vector<vector<vector<double> > >& spot_intensities_,
	             const vector<Vector<4> >& spots_,
	             const vector<vector<int> >& spot_pixels_,
	             const Matrix<3> A_,
	             const Vector<3> pi_,
	             double variance_,
	             int sample_iterations_)
	:pixel_intensities(pixel_intensities_), //pixel_intensities: [frame][pixels]
	 spot_intensities(spot_intensities_),   //spot_intensities: [spot][frame][pixel]
	 spots(spots_),
	 spot_pixels(spot_pixels_),
	 A(A_),
	 pi(pi_),
	 base_variance(variance_),
	 variance(variance_),
	 sample_iterations(sample_iterations_),
	 num_frames(pixel_intensities.size()),
	 num_pixels(pixel_intensities[0].size()),
		//Observations vector. As usual for this application, the observations are just
		//numbered integers which refer to data held elsewhere.
	 O(sequence(num_frames)),
		//Start all spots OFF, so the intensity is 0. OFF is 1 or 2, not 0!!!
		//sample_list: [sample][spot][frame]: list of samples drawn using Gibbs sampling
	 current_sample(spots.size(), vector<State>(num_frames, 2)), //current sample [spot][frame]
		//pixel intensities assosciated with the current sample [frame][pixel]
	 current_sample_intensities(num_frames, vector<double>(num_pixels))
	{
		//Check a bunch of stuff
		assert_same_size(pixel_intensities);
		assert_same_size(spot_intensities);
	}
	
	///Update the noide variance. Used for adding thermal noise.
	///@param v noise variance.
	void set_variance(double v)
	{
		variance = v;
	}

	///Reset the gibbs sampler oro the initial state (all spots off)
	void reset()
	{
		vector<State> off(num_frames, 2);
		fill(current_sample.begin(), current_sample.end(), off);

		vector<double> black(num_pixels);
		fill(current_sample_intensities.begin(), current_sample_intensities.end(), black);
		variance = base_variance;
	}
	
	///Get the next sample
	///@param rng Random number generator
	template<class T> void next(T& rng)
	{
	
		std::vector<array<double, 3> > delta3;
		for(int j=0; j < sample_iterations; j++)
			for(int k=0; k < (int) spots.size(); k++)
			{
				//Subtract off the spot we're interested in.
				remove_spot(current_sample_intensities, spot_intensities[k], current_sample[k], spot_pixels[k]);

				//Now current_sample_intensities is the image value for every spot in every frame,
				//except the current spot, which is always set to off. This allows us to add it in 
				//easily.

				SpotWithBackgroundDriftMasked B3(current_sample_intensities, spot_intensities[k], pixel_intensities, variance, spot_pixels[k]);
				forward_algorithm_delta2<3>(A, pi, B3, O, delta3);

				current_sample[k] = backward_sampling<3,State, T>(A, delta3, rng);
				//Put the newly sampled spot in
				add_spot(current_sample_intensities, spot_intensities[k], current_sample[k], spot_pixels[k]);
			}
	}
/*	void next()
	{
		RngDrand48 rng;
		next(rng);
	}	
*/
	///Retrieve the current sample
	const vector<vector<State> >& sample() const
	{
		return current_sample;
	}
	///Retrieve the intensities for the current sample
	const vector<vector<double> >& sample_intensities() const
	{
		return current_sample_intensities;
	}

};
#endif 

}

using SampledMultispot::SpotWithBackground;
using SampledMultispot::remove_spot;
using SampledMultispot::add_spot;
using SampledMultispot::compute_spot_intensity;
using SampledMultispot::compute_spot_intensity_hessian;
using SampledMultispot::compute_spot_intensity_derivatives;
using SampledMultispot::sequence;
using SampledMultispot::GibbsSampler;
using SampledMultispot::GibbsSampler2;

#endif
