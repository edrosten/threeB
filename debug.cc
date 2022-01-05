#include <cvd/byte.h>

/** Scales an image in to the correct range for bytes.
@param hi Brightest pixel in the image 
@param lo Dimmest pixel in the image
@param im Image to scale
@returns scaled image
@ingroup gDebug
*/
Image<byte> scale_to_bytes(const Image<float>& im, float lo, float hi)
{
	Image<byte> out(im.size());
	for(int r=0; r < out.size().y-0; r++)
		for(int c=0; c < out.size().x-0; c++)
			out[r][c] = (int)floor((im[r][c]-lo)*255/(hi-lo));

	return out;
}

/** Find the variance of every patch in the image and save it to a file
@ingroup gDebug
@param ims List of images.
*/
void test_output_patch_variance(const vector<Image<float> >&  ims)
{
	assert_same_size(ims);

	int rr = GV3::get<int>("test.variance.radius", 1, -1);
	ImageRef r(rr, rr);
	ImageRef size = r*2 + ImageRef(1,1);

	Image<float> stds(ims.front().size(), 0);

	ImageRef p;	
	for(ImageRef p(0,0); p.y < stds.size().y - size.y; p.y++)
	{
		for(p.x=0; p.x < stds.size().x - size.x; p.x++)
			stds[p + r] = sqrt(mean_and_variance(sub_images(ims, p, size)).second);
	}

	SubImage<float> s = stds.sub_image(ImageRef(2,2), stds.size() - ImageRef(4,4));

	float hi = *max_element(s.begin(), s.end());
	float lo = *min_element(s.begin(), s.end());
	cerr << hi << " " << lo << endl;
	img_save(scale_to_bytes(stds, lo, hi), "test_variance.png");
}





