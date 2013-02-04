#ifndef INC_DRIFT_H
#define INC_DRIFT_H


///Compute the integer step number from the frame number
///
///Drift is linear, with a piecewise constant approximation
///@ingroup gMultiSpotDrift
///@param frames Number of frames
///@param steps Number of steps
///@param frame Current frame number
int frame_to_step(int frames, int steps, int frame)
{
	return (frame * steps) / frames;
}


///Compute the approximate, real-valued frame number from the step number
///
///Drift is linear, with a piecewise constant approximation
///@ingroup gMultiSpotDrift
///@param frames Number of frames
///@param steps Number of steps
///@param frame Current frame number
double step_to_approximate_frame(int frames, int steps, int step)
{
	return (step + 0.5)/steps * frames;
}


///Apply the drift model to the spot
///@ingroup gMultiSpotDrift
///@param spot Spot parameters (size/brightness/position)
///@param drift Drift vector
///@param frame_number Frame number, note that it is real-valued.
Vector<4> drift_spot(const Vector<4>& spot, const Vector<2>& drift, double frame_number)
{
	return spot + frame_number * makeVector(0, 0, drift[0], drift[1]);
}


#endif
