#include <iomanip>

#include <tag/printf.h>
#undef make_tuple
#include <gvars3/instances.h>

#include <cvd/image_io.h>
#include <cvd/draw.h>
#include "forward_algorithm.h"
#include "storm.h"
#include "mt19937.h"

using namespace std;
using namespace CVD;
using namespace GVars3;
using namespace TooN;
using namespace std::tr1;

template<int States> vector<int> run_hmm(Matrix<States> A, Vector<States> pi, int n, MT19937& rng)
{

	bool bad=1;

	vector<int> states;

	while(bad)
	{
		states.clear();
	
		int state = select_random_element(pi, 1.0, rng);

		for(int i=0; i < n;i++, state = select_random_element(A[state], 1.0, rng))
		{
			states.push_back(state);

			if(state==0)
				bad=0;
		}
	}

	return states;
}

double find_median_threshold(const vector<ImageRef>& disc, const Image<float>& im)
{
	//Find the threshold at which a median filter would detect a spot
	vector<float> pix(disc.size());

	for(unsigned int i=0; i < disc.size(); i++)
		pix[i] = im[disc[i]];
	
	nth_element(pix.begin(), pix.begin() + pix.size()/2, pix.end());

	return *(pix.begin() + pix.size()/2);
}


int main(int argc, char ** argv)
{
	GUI.LoadFile("simulations.cfg");
	GUI.parseArguments(argc, argv);

	ImageRef size = GV3::get<ImageRef>("simulations.size", "", -1);
	int num_frames = GV3::get<int>("simulations.frames", 0, -1);
	int num_spots = GV3::get<int>("simulations.spots", 0, -1);
	double spot_max_brightness = GV3::get<double>("simulations.max_brightness", 0, -1);
	double spot_fwhm = GV3::get<double>("simulations.fwhm", 0, -1);
	double spot_sigma = spot_fwhm / ( 2 * sqrt(2 * log(2)));
	double spot_brightness = spot_max_brightness * spot_sigma * sqrt(2*M_PI);

	double spot_spacing = GV3::get<double>("simulations.spacing", 0, -1);
	
	int seed = GV3::get<int>("simulations.seed", 0 , 1);
	bool find_threshold=GV3::get<bool>("simulations.find_threshold", 0, 1);

	//This really only makes sense for a single spot
	if(num_spots != 1)
		find_threshold=0;
	
	vector<vector<ImageRef> > discs;

	if(find_threshold)
	{
		vector<float> disc_radii = GV3::get<vector<float> >("simulations.discs", "", -1);
		for(unsigned int i=0; i < disc_radii.size(); i++)
		{
			discs.push_back(getDisc(disc_radii[i]));
			for(unsigned int i=0; i < discs.back().size(); i++)
				discs.back()[i] += size/2;
		}
	}

	MT19937 rng;
	
	rng.simple_seed(seed);

	string stub = GV3::get<string>("stub", "img_%06i.tif");
	string stub2 = GV3::get<string>("stub2", "");

	Matrix<3> A  = GV3::get<Matrix<3> >("A", Zeros, 1);
	Vector<3> pi = GV3::get<Vector<3> >("pi", Zeros, 1);
	

	//Generate the hidden states
	vector<vector<int> >  states; // [spot][frame]
	for(int s=0; s < num_spots; s++)
		states.push_back(run_hmm(A, pi, num_frames, rng));

	vector<Vector<4> > spot_parameters;
	
	cout << "[ ";
	for(int s=0; s < num_spots; s++)
	{
		double xpos = size.x / 2.0 + (s - (num_spots-1)/2.0) * spot_spacing;
		spot_parameters.push_back(makeVector(spot_brightness, spot_sigma, xpos, size.y/2.0));
		cout << spot_parameters.back() << " ";
	}
	cout << " ]\n";

	Image<float> out(size);
	

	for(int f=0; f < num_frames; f++)
	{
		for(int y=0; y < size.y; y++)
			for(int x=0; x < size.x; x++)
			{
				out[y][x] = rng.gaussian();

				for(int s=0; s < num_spots; s++)
					if(states[s][f] == 0)
						out[y][x] += spot_shape(makeVector(x,y), spot_parameters[s]);
			}

		img_save(out, sPrintf(stub, f));

		if(find_threshold)
		{
			cout << "TH " << states[0][f];
			for(unsigned int i=0; i < discs.size(); i++)
				cout << " " << scientific << find_median_threshold(discs[i], out);
			cout << endl;
		}

		if(stub2 != "")
		{
			for(int y=0; y < size.y; y++)
				for(int x=0; x < size.x; x++)
					out[y][x] = max(0.,min(1., out[y][x] / 10 + .2));

			img_save(out, sPrintf(stub2, f));
		}

	}
}
