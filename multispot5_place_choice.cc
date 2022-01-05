#include "multispot5_place_methods.h"
#include "multispot5_place_choice.h"
#include "debug.h"
#include <gvars3/instances.h>
#include <algorithm>
#include <cvd/morphology.h>
#include <cvd/draw.h>
#include <cvd/byte.h>

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
		float dilate_mask_radius = GV3::get<float>("dilate_mask_as_filter_radius", 0.0, 1);
		vector<ImageRef> mask;


		if(dilate_mask_radius != 0)
		{
			Image<byte> filter(ims[0].size());
			filter.fill(0);
			for(unsigned int i=0; i < region.size(); i++)
				filter[region[i]]=1;

			Image<byte> dilated = morphology(filter, getDisc(dilate_mask_radius), Morphology::BinaryDilate<byte>());

			for(int r=0; r < dilated.size().y; r++)
				for(int c=0; c < dilated.size().x; c++)
					if(dilated[r][c])
						mask.push_back(ImageRef(c, r));
			
			GV3::get<vector<ImageRef> >("filter", "", 1) = region;
		}
		else
		{
			mask=region;
		}


		string placement = GV3::get<string>("placement", "uniform", 1);

	
		save_spots << extra << endl;

		if(placement == "ye_olde")
		{
			StateParameters p(generate_state_parameters_ye_olde(log_ratios, ims, mask));
			fit_spots_new(ims, p, save_spots, g, ui);
		}
		else if(placement == "uniform")
		{
			int num = GV3::get<int>("placement.uniform.num_spots", 0, -1);
			StateParameters p(place_spots_uniform(num, mask, log_ratios.size()));
			fit_spots_new(ims, p, save_spots, g, ui);
		}
		else if(placement == "intensity_sampled")
		{
			int num = GV3::get<int>("placement.uniform.num_spots", 0, -1);
			StateParameters p(place_spots_intensity_sampled(num, mask, ims));
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

