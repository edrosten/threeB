#include <tag/printf.h>
#undef make_tuple

#include <tr1/tuple>
#include <algorithm>
#include <climits>
#include <iomanip>
#include <map>
#include <cvd/image_io.h>
#include <cvd/image_convert.h>
#include <cvd/morphology.h>
#include <cvd/connected_components.h>
#include <cvd/draw.h>
#include <cvd/vector_image_ref.h>
#include <cvd/byte.h>

#include <gvars3/instances.h>

#include "storm_imagery.h"
#include "multispot5.h"
#include "multispot5_place_choice.h"
#include "utility.h"

using namespace std;
using namespace std::tr1;
using namespace CVD;
using namespace GVars3;
using namespace TooN;



vector<vector<ImageRef> > get_regions(const SubImage<double>& log_ratios)
{
	gvar3<double> radius("radius", 0, 1);

	//Set the liklihood ratio threshold/spot density prior
	//same thing.
	double threshold = GV3::get<double>("threshold", 0, -1);


	//Threshold image
	Image<byte> thresholded(log_ratios.size(), 0);
	for(int r=0; r < thresholded.size().y; r++)
		for(int c=0; c < thresholded.size().x; c++)
			thresholded[r][c] = 255 * (log_ratios[r][c] > threshold);
	
	//Dilate
	Image<byte> dilated = morphology(thresholded, getDisc(*radius), Morphology::BinaryDilate<byte>());

	transform(dilated.begin(), dilated.end(), dilated.begin(), bind1st(multiplies<int>(), 255));
	
	//Connected components of dilated image
	vector<ImageRef> fg;
	for(int r=0; r < thresholded.size().y; r++)
		for(int c=0; c < thresholded.size().x; c++)
			if(dilated[r][c])
				fg.push_back(ImageRef(c, r));

	vector<vector<ImageRef> > regions;
	connected_components(fg, regions);

	return regions;
}

void mmain(int argc, char** argv)
{
	GUI.LoadFile("multispot5.cfg");
	int lastarg = GUI.parseArguments(argc, argv);
	if(lastarg >= argc)
	{	
		cerr << "Specify the images to load\n";
		exit(1);
	}
	vector<string> files(argv + lastarg, argv + argc);
	
	//Save this now since the de-checkpointing code will kl0bber it 
	//when it reloads the gvars
	string save_spots_file = GV3::get<string>("save_spots", "", -1);

	string checkpoint_file = GV3::get<string>("load_checkpoint", "", 1);

	if(checkpoint_file != "")
	{
		//Load and de-checkpointing
		ifstream chk;
		open_or_die(chk, checkpoint_file);
		
		StateParameters p;

		try{
			p = parse_log_file(chk);
		}
		catch(LogFileParseError e)
		{
			cerr << "SI TEH FUX0R11ONEone!oneleven: " << e.what << endl;
			exit(1);
		}
		
		vector<Image<float> > ims = load_and_normalize_images(files);

		//Restore kl0bbered variable
		GV3::get<string>("save_spots") = save_spots_file;
		
		ofstream save_spots;
		open_or_die(save_spots, save_spots_file);
		
		fit_spots_new(ims, p, save_spots, *null_graphics());

	}
	else
	{
		vector<Image<float> > ims = load_and_normalize_images(files);

		//Load the log_ratios image.
		//We will use this as a starting point for searching for spots.
		Image<double> log_ratios;
		try
		{
			log_ratios = img_load(GV3::get<string>("log_ratios", "", -1));
		}
		catch(Exceptions::All e)
		{
			cerr << "Error loading " << GV3::get<string>("log_ratios", "") << ": " << e.what << endl;
			exit(1);
		}


		gvar3<int> cluster_to_show("cluster_to_show", 0, -1);
		gvar3<int> use_largest("use_largest", 0, 1);

		vector<vector<ImageRef> > regions;

		regions = get_regions(log_ratios);
		if(regions.size() == 0)
		{
			cerr << "There are no regions!\n";

			ofstream save_spots;
			open_or_die(save_spots, save_spots_file);
			save_spots << "NOREGIONS\n";

			exit(1);
		}
		
		if(*use_largest && !regions.empty())
		{
			*cluster_to_show=0;
			for(unsigned int i=1; i < regions.size(); i++)
				if(regions[i].size() > regions[*cluster_to_show].size())
					*cluster_to_show = i;
					
		}
		else
			*cluster_to_show = max(min(*cluster_to_show, (int)regions.size() - 1), 0);

		
		auto_ptr<FitSpotsGraphics> gr = null_graphics();
		place_and_fit_spots(ims, regions[*cluster_to_show], log_ratios, save_spots_file, *gr);
	}
}
	
int main(int argc, char** argv)
{
	try{
		mmain(argc, argv);
	}
	catch(Exceptions::All e)
	{
		cerr << "Fatal error: " << e.what << endl;
	}
}
