#include "multispot5_place_methods.h"
#include "multispot5_place_choice.h"
#include "debug.h"
#include <gvars3/instances.h>
#include <algorithm>

using namespace std;
using namespace CVD;
using namespace GVars3;
using namespace TooN;

void place_and_fit_spots(const vector<Image<float> >& ims, const vector<ImageRef>& region_, const Image<double>& log_ratios, string save_spots_file, FitSpotsGraphics& g, const string& extra)
{
	auto_ptr<UserInterfaceCallback> ui(null_ui());
	ofstream save_spots;
	open_or_die(save_spots, save_spots_file);
	place_and_fit_spots(ims, region_, log_ratios, save_spots, g, *ui, extra);
}

void place_and_fit_spots(const vector<Image<float> >& ims, const vector<ImageRef>& region_, const Image<double>& log_ratios, ofstream& save_spots, FitSpotsGraphics& g, UserInterfaceCallback& ui, const std::string& extra)
{	
	assert_same_size(ims);
	assert(ims[0].size() == log_ratios.size());

	vector<ImageRef> region = region_;
	sort(region.begin(), region.end());

	string mode = GV3::get<string>("mode", "new", 1);
	
	if(mode == "new")
	{
		string placement = GV3::get<string>("placement", "uniform", 1);

	
		save_spots << extra << endl;

		if(placement == "ye_olde")
		{
			StateParameters p(generate_state_parameters_ye_olde(log_ratios, ims, region));
			fit_spots_new(ims, p, save_spots, g, ui);
		}
		else if(placement == "uniform")
		{
			int num = GV3::get<int>("placement.uniform.num_spots", 0, -1);
			StateParameters p(place_spots_uniform(num, region, log_ratios.size()));
			fit_spots_new(ims, p, save_spots, g, ui);
		}
		else if(placement == "intensity_sampled")
		{
			int num = GV3::get<int>("placement.uniform.num_spots", 0, -1);
			StateParameters p(place_spots_intensity_sampled(num, region, ims));
			fit_spots_new(ims, p, save_spots, g, ui);
		}
		else
			cerr << "Mode must be uniform or ye_olde, not `" + placement + "'.\n";
	}
	else
	{
		cerr << "Mode must be new , not `" + mode + "'.\n";
	}

}

