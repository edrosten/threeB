#include "utility.h"
#include "debug.h"
#include <climits>
using namespace std;
using namespace CVD;

//! @cond Doxygen_Suppress
const std::vector<CVD::SubImage<float> > sub_images(const std::vector<CVD::Image<float> >& im, CVD::ImageRef pos, ImageRef size)
{
	assert_same_size(im);

	vector<SubImage<float> > subs;

	for(unsigned int i=0; i < im.size(); i++)
		subs.push_back(im[i].sub_image(pos, size));
	return subs;
}

pair<ImageRef, ImageRef> boundingbox(const vector<ImageRef> & all_spots)
{
	ImageRef lo(INT_MAX, INT_MAX), hi(INT_MIN, INT_MIN);
	for(unsigned int i=0; i < all_spots.size(); i++)
	{
		lo[0] = min(lo[0], all_spots[i][0]);
		lo[1] = min(lo[1], all_spots[i][1]);

		hi[0] = max(hi[0], all_spots[i][0]);
		hi[1] = max(hi[1], all_spots[i][1]);
	}

	return make_pair(lo, hi - lo + ImageRef(1,1));
}

//! @endcond 
