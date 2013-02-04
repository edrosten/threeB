#ifndef MULTISPOT5_H
#define MULTISPOT5_H
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <tr1/memory>
#include <cvd/image.h>
#include <cvd/byte.h>
#include <TooN/TooN.h>
#include <TooN/so2.h>

#include "utility.h"
#include "mt19937.h"


/// Graphics class for FittingSpots.
/// This abstraction prevents FitSpots from depending on and graphics library.
/// The functions are tied very heavily to the internals of FitSpots.
/// @ingroup gStorm
class FitSpotsGraphics
{
	public:
		///Initialize graphics
		///@param size Size of window for display
		virtual void init(CVD::ImageRef size)=0;

		///Draw a bunch of stuff
		///@param spots List of spots to draw
		///@param im Background image
		///@param box Bounding box of region
		///@param N Spot to highlight
		///@param s Extra spot to draw as a cross
		virtual void draw_krap(const std::vector<TooN::Vector<4> >& spots, const CVD::Image<CVD::byte>& im, const BBox& box, int N, TooN::Vector<4> s = TooN::Ones * 1e99)=0;

		///Swap buffers if double buffered
		virtual void swap()=0;

		///Draw the pixel mask in an (r,g,b,a) tuple colour
		///@param pix mask
		///@param r red
		///@param g green
		///@param b blue
		///@param a alpha
		virtual void draw_pixels(const std::vector<CVD::ImageRef>& pix, float r, float g, float b, float a=1)=0;

		///Draw a bounding box
		///@param bbox Box corners
		virtual void draw_bbox(const BBox& bbox)=0;

		///Draw a cross
		///@param p Position of cross
		///@param size Size to draw cross
		virtual void glDrawCross(const TooN::Vector<2>& p, int size=3)=0;

		///Desctructor
		virtual ~FitSpotsGraphics();
};


/// Callback class used by FitSpots to provide enough hooks for a user interface.
/// @ingroup gStorm
class UserInterfaceCallback
{
	public:
		///This function is called once per spot in each pass. The idea is to
		///provide a display along the lines of: Iteration #1 optimizing #2% complete
		///@param iteration Iteration number
		///@param pass      Pass number
		///@param spot_num  Spot currently being optimized
		///@param total_spots Total number of spots to be optimized
		virtual void per_spot(int iteration, int pass, int spot_num, int total_spots)=0;

		///This function is called once per spot in the modification phase. The idea is to
		///provide a display along the lines of: Iteration #1 modifying #2% complete
		///@param iteration Iteration number
		///@param spot_num  Spot currently being optimized
		///@param total_spots Total number of spots to be optimized
		virtual void per_modification(int iteration, int spot_num, int total_spots)=0;

		///This function is called once each time PASS data is outputted.
		///It will allow the GUI to build up a reconstruction.
		///@param iteration Iteration number
		///@param pass      Pass number
		///@param spots     Data to be reconstructed
		virtual void per_pass(int iteration, int pass, const std::vector<TooN::Vector<4> >& spots)=0;

		///The user wishes to issue a stop instruction to the program (perhaps done via an 
		///asynchronus call to an instance of of UserInterfaceCallback). This function is 
		///called as often as possible and will throw UserIssuedStop when the condition is
		///met.
		virtual void perhaps_stop()=0;
		
		virtual ~UserInterfaceCallback();
		struct UserIssuedStop{};
};

std::auto_ptr<FitSpotsGraphics> null_graphics();
std::auto_ptr<UserInterfaceCallback> null_ui();

///Null struct thrown if a parse error is encountered when trying to load a log file.
struct LogFileParseError
{
	///@param s error message to set
	LogFileParseError(const std::string &s)
	:what(s)
	{}
	
	/// variable holding the parse error error message
	std::string what;
};

///Internal state (excluding fixed settings) which represents the entire
///internal state of spot fitting. Used to restart from interruptions.
struct StateParameters{
	std::tr1::shared_ptr<MT19937> rng;    ///< Random number generator state
	std::vector<TooN::Vector<4> > spots;  ///< Spots positions 
	int pass;                             ///< Pass number
	int iteration;                        ///< Iteration number
	std::vector<CVD::ImageRef> pixels;    ///< Area for analysis
};

StateParameters generate_state_parameters_ye_olde(const CVD::BasicImage<double>& log_ratios, const std::vector<CVD::Image<float> >& ims, std::vector<CVD::ImageRef> pixels);

void fit_spots_new(const std::vector<CVD::Image<float> >& ims, StateParameters& p, std::ofstream& save_spots, FitSpotsGraphics&);
void fit_spots_new(const std::vector<CVD::Image<float> >& ims, StateParameters& p, std::ofstream& save_spots, FitSpotsGraphics&, UserInterfaceCallback&);

StateParameters parse_log_file(std::istream& in);
#endif
