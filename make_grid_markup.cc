#include <cvd/image_io.h>
#include <cvd/draw.h>
#include <cvd/morphology.h>
#include <gvars3/instances.h>
#include <tag/printf.h>
#include <climits>

using namespace CVD;
using namespace std;
using namespace tag;
using namespace GVars3;

int main(int argc, char** argv)
{	
	GUI.parseArguments(argc, argv);

	float dilate = 0;
	int size = GV3::get<int>("size", 0, -1);

	Image<byte> im = img_load(GV3::get<string>("image"));

	int d = ceil(dilate);

	vector<ImageRef> disc = getDisc(dilate);

	Image<byte> mask(im.size()), filter(im.size()), final(im.size());

	final.zero();
	int n=0;

	string mask_str = GV3::get<string>("mask", "", 1);

	if(mask_str != "")
	{
		dilate = GV3::get<float>("radius", 0, -1);
	}

	string filt_str = GV3::get<string>("filter", "", -1);
	
	int best_cells_filled = INT_MAX, bestxoff=0, bestyoff = 0;
	
	//Find the optiml offset
	for(int yoff=0; yoff < size; yoff++)
		for(int xoff=0; xoff < size; xoff++)
		{
			int cells_filled=0;
			for(int y=yoff; y < im.size().y-size; y += size)
				for(int x=xoff; x < im.size().x-size; x += size)
				{
					if(x-d >= 0 && x+d < im.size().x && y-d >= 0 && y+d < im.size().y)
					{
						for(int i=0; i < size; i++)
							for(int j=0; j < size; j++)
							{
								if(im[y+j][x+i])
								{
									cells_filled++;
									goto cont;
								}
							}
					}	
					cont:;
				}

			clog << "Testing " << xoff << ", " << yoff << " " << cells_filled << "\n";

			if(cells_filled < best_cells_filled)
			{
				best_cells_filled = cells_filled;
				bestxoff = xoff;
				bestyoff = yoff;
			}
		}

	clog << "---------\n";
	clog << "Using " << best_cells_filled << ", " << bestyoff << " --> " << best_cells_filled << endl;

	for(int y=bestyoff; y < im.size().y-size; y += size)
		for(int x=bestxoff; x < im.size().x-size; x += size)
			if(x-d >= 0 && x+d < im.size().x && y-d >= 0 && y+d < im.size().y)
			{
				filter.zero();
				mask.zero();

				bool use=0;

				for(int i=0; i < size; i++)
					for(int j=0; j < size; j++)
					{
						if(im[y+j][x+i])
							use=1;
						filter[y+j][x+i] = 255;
					}
				
				if(use)
				{
					for(int i=0; i < size; i++)
						for(int j=0; j < size; j++)
							final[y+j][x+i] = 255;

					img_save(filter, sPrintf(filt_str, n));

					if(mask_str != "")
					{
						morphology(filter, disc, Morphology::Dilate<byte>(), mask);
						img_save(mask, sPrintf(mask_str, n));
					}

					n++;

				}
			}


	img_save(final, "final_xxxxx.png");
}
