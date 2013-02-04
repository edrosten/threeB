#include <gvars3/instances.h>
#include <cvd/image_io.h>
#include <cvd/convolution.h>
#include <TooN/wls.h>
#include <tr1/tuple>

#include "storm_imagery.h"
#include "debug.h"
#include "utility.h"

using namespace CVD;
using namespace TooN;
using namespace GVars3;
using namespace std;
using namespace std::tr1;

/**Load all images from disk and do the initial preprocessing. 

@param names List of filenames to load.
@return preprocessed images.
@ingroup gStormImages
**/
vector<Image<float> > load_and_preprocess_images2(const vector<string>& names)
{
	vector<Image<float> > ims;
	//Load images
	for(unsigned int i=0; i < names.size(); i++)
	{
		Image<float> im = img_load(names[i]);
		ims.push_back(im);

		if(ims.back().size() != ims[0].size())
		{
			cerr << "Error with image " << names[i] << ":  all images must be the same size!\n";
			exit(1);
		}
	}
	double mean, variance;
	tie(mean, variance) = mean_and_variance(ims);
	{
		for(unsigned int i=0; i < ims.size(); i++)
			transform(ims[i].begin(), ims[i].end(), ims[i].begin(), bind2nd(minus<double>(), mean));
		for(unsigned int i=0; i < ims.size(); i++)
			transform(ims[i].begin(), ims[i].end(), ims[i].begin(), bind1st(multiplies<double>(), 1/ sqrt(variance)));
	}
	
	tie(mean, variance) = mean_and_variance(ims);

	cerr << "Rescaled:\n";
	cerr << "mean = " << mean << endl;
	cerr << "std  = " << sqrt(variance) << endl;



	//Normalize...

	//Fit the background model
	ImageRef size = ims[0].size();
	Vector<10> p = Zeros;	
	p[6]=-3;
	p[9]=-4;

	Image<Vector<6> > monomials(size);
	Image<double> polynomial(size);
	for(int yy=0; yy < size.y; yy++)
		for(int xx=0; xx < size.x; xx++)
		{
			double x = xx *2./ size.x -1 ;
			double x2 = x*x;
			double y = yy *2./size.y - 1;
			double y2 = yy;
			monomials[yy][xx] = makeVector(1, x, y, x2, x*y, y2);
		}
					

	for(int i=0;i < 100;i++)
	{
		for(int yy=0; yy < size.y; yy++)
			for(int xx=0; xx < size.x; xx++)
				polynomial[yy][xx] = monomials[yy][xx] *  p.slice<0,6>();

		WLS<10> wls;	
		for(unsigned int i=0; i < ims.size(); i++)
			for(int yy=0; yy < size.y; yy++)
				for(int xx=0; xx < size.x; xx++)
				{
					double t = i *1. / ims.size();
					double func = polynomial[yy][xx] * (exp(p[6]*t) + p[8]*exp(p[9]*t)) + p[7];

					Vector<10> mJ;

					mJ.slice<0,6>() = exp(p[6]*t)* monomials[yy][xx];
					//mJ.slice<3,3>() = Zeros;
					mJ[6] = polynomial[yy][xx] * exp(p[6]*t) * t;
					//mJ[6] = func  * t;
					mJ[7] = 1;

					mJ[8] = polynomial[yy][xx] * exp(p[9]*t);
					mJ[9] = polynomial[yy][xx] * exp(p[9]*t) * t * p[8];

					double err = ims[i][yy][xx] - func;

					double w;

					
					if(err > 0)
						w = .01 / (abs(err) + .01);
					else
						w = 1;

					wls.add_mJ(func - ims[i][yy][xx], -mJ, w);
				}
		
		wls.add_prior(10);
		wls.compute();

		p += wls.get_mu();

		cout << p << endl << endl;
	}
	
	for(unsigned int i=0; i < ims.size(); i++)
		for(int yy=0; yy < size.y; yy++)
			for(int xx=0; xx < size.x; xx++)
			{
				double x = xx *2./ size.x -1 ;
				double x2 = x*x;
				double y = yy *2./size.y - 1;
				double y2 = yy;
				double t = i *1. / ims.size();
				Vector<6> f = makeVector(1, x, y, x2, x*y, y2);
				
				double func = f * p.slice<0,6>() * (exp(p[6]*t) + p[8]*exp(p[9]*t)) + p[7];
				ims[i][yy][xx] -= func;
			}

	tie(mean, variance) = mean_and_variance(ims);

	//A sanity check.
	cerr << "The mean should be small compared to std:\n";
	cerr << "mean = " << mean << endl;
	cerr << "std  = " << sqrt(variance) << endl;

	//Scale by the variance.
	{
		for(unsigned int i=0; i < ims.size(); i++)
			transform(ims[i].begin(), ims[i].end(), ims[i].begin(), bind1st(multiplies<double>(), 1/ sqrt(variance)));
	}
	tie(mean, variance) = mean_and_variance(ims);

	cerr << "Rescaled:\n";
	cerr << "mean = " << mean << endl;
	cerr << "std  = " << sqrt(variance) << endl;

	return ims;
}




/**Load all images from disk and do the initial preprocessing. Currently
this is a high pass filter to make the resultimg images zero mean.

The filter is controlled with the \c preprocess.lpf and \c preprocess.skip Gvars

See also load_and_preprocess_image()

@param names List of filenames to load.
@return preprocessed images.
@ingroup gStormImages
**/
vector<Image<float> > load_and_preprocess_images(const vector<string>& names)
{
	vector<Image<float> > ims;

	//float wide = GV3::get<float>("preprocess.lpf", 0., -1);
	//bool p = GV3::get<bool>("preprocess.skip", 0, -1);
	
	for(unsigned int i=0; i < names.size(); i++)
	{
		Image<float> im = img_load(names[i]);
	
		ims.push_back(preprocess_image(im));
	
		if(ims.back().size() != ims[0].size())
		{
			cerr << "Error with image " << names[i] << ":  all images must be the same size!\n";
			exit(1);
		}
	}
	return ims;
}

/**Compute the mean and variance of the (on average) darkest pixels, in order
to find the correct scaling, by examining hte background.
*/
pair<double, double> auto_fixed_scaling(const vector<Image<float> >& ims, double frac)
{
	assert_same_size(ims);
	
	//Compute the mean image (ish)
	Image<double> ave(ims[0].size());
	ave.fill(0);
	for(unsigned int i=0; i < ims.size(); i++)
		for(int y=0; y < ave.size().y; y++)
			for(int x=0; x < ave.size().x; x++)
				ave[y][x] += ims[i][y][x];
	
	//Find the smallest N% of the pixels...
	vector<pair<double, ImageRef> > pixels;
	for(int y=0; y < ave.size().y; y++)
		for(int x=0; x < ave.size().x; x++)
			pixels.push_back(make_pair(ave[y][x], ImageRef(x,y)));

	int npix = (int) floor(frac *pixels.size() + 0.5);
	npix = max(0, min(npix, (int) pixels.size()));

	nth_element(pixels.begin(), pixels.begin() + npix, pixels.end());

	pixels.resize(npix);
	
	//Now compute the mean and variance of those pixels.
	double sum=0, sum2=0;

	for(unsigned int i=0; i < ims.size(); i++)
	{	
		for(unsigned int j=0; j < pixels.size(); j++)
		{
			sum += ims[i][pixels[j].second];
			sum2 += sq(ims[i][pixels[j].second]);
		}
	}

	double num = 1.0 * pixels.size() * ims.size();
	double mean = sum / num;
	double std  = sqrt(((sum2/num) - mean*mean) * num / (num-1));

	cout << "Automatic determination of fixed scaling:" << endl
	     << "mean       = " << mean << endl
		 << "std        = " << std << endl
		 << "sqrt(mean) = " << sqrt(mean*255)/255 << endl;
	
	return make_pair(mean, std);
}

/**Wrapper for load_and_preprocess_images() to allow more flexible behaviour.

@param files List of filenames to load.
@return preprocessed images.
@ingroup gStormImages
**/
vector<Image<float> > load_and_normalize_images(const vector<string>& files)
{
	//Load the raw data, and then load the spot parameters.
	vector<Image<float> > ims = load_and_preprocess_images(files);
	double mean, variance;
	tie(mean, variance) = mean_and_variance(ims);

	if(GV3::get<bool>("preprocess.fixed_scaling", 0, FATAL_IF_NOT_DEFINED))
	{
		bool skip = GV3::get<bool>("preprocess.skip");
		if(!skip)
		{
			cerr << "WARNING WARNING WARNING WARNING!!!!!!!!!!!!!!!\n";
			cerr << "preprocessing and fixed scaling selected!!!\n";
			exit(1);
		}

		double sub, div;
		if(GV3::get<bool>("preprocess.fixed_scaling.auto", 0, FATAL_IF_NOT_DEFINED))
		{
			double frac = GV3::get<double>("preprocess.fixed_scaling.auto.proportion", 0, FATAL_IF_NOT_DEFINED);
			tie(sub, div) = auto_fixed_scaling(ims, frac);
		}
		else
		{
			sub = GV3::get<double>("preprocess.fixed_scaling.subtract", 0, FATAL_IF_NOT_DEFINED);
			div = GV3::get<double>("preprocess.fixed_scaling.divide", 0, FATAL_IF_NOT_DEFINED);
		}

		for(unsigned int i=0; i < ims.size(); i++)
			for(Image<float>::iterator j=ims[i].begin(); j != ims[i].end(); j++)
				*j = (*j - sub)/div;
	}
	else
	{
		//A sanity check.
		cerr << "The mean should be small compared to std:\n";
		cerr << "mean = " << mean << endl;
		cerr << "std  = " << sqrt(variance) << endl;

		//Scale by the variance.
		{
			for(unsigned int i=0; i < ims.size(); i++)
				transform(ims[i].begin(), ims[i].end(), ims[i].begin(), bind1st(multiplies<double>(), 1/ sqrt(variance)));
		}
	}

	tie(mean, variance) = mean_and_variance(ims);

	//A sanity check.
	cerr << "Rescaled:\n";
	cerr << "mean = " << mean << endl;
	cerr << "std  = " << sqrt(variance) << endl;

	return ims;
}



