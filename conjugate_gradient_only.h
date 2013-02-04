#ifndef CONJUGATE_GRADIENT_ONLY
#define CONJUGATE_GRADIENT_ONLY
#include <TooN/TooN.h>
#include <utility>
#include <cstdlib>

///Class for performing optimization with Conjugate Gradient, where only the derivatives are available.
///@ingroup gStorm
template<int Size=-1, class Precision=double> struct ConjugateGradientOnly
{
	const int size;      ///< Dimensionality of the space.
	TooN::Vector<Size> g;      ///< Gradient vector used by the next call to iterate()
	TooN::Vector<Size> new_g;  ///< The new gradient set at the end of iterate()
	TooN::Vector<Size> h;      ///< Conjugate vector to be searched along in the next call to iterate()
	TooN::Vector<Size> minus_h;///< negative of h as this is required to be passed into a function which uses references (so can't be temporary)
	TooN::Vector<Size> old_g;  ///< Gradient vector used to compute $h$ in the last call to iterate()
	TooN::Vector<Size> old_h;  ///< Conjugate vector searched along in the last call to iterate()
	TooN::Vector<Size> x;      ///< Current position (best known point)
	TooN::Vector<Size> old_x;  ///< Previous point (not set at construction)

	Precision tolerance; ///< Tolerance used to determine if the optimization is complete. Defaults to square root of machine precision.
	int       max_iterations; ///< Maximum number of iterations. Defaults to \c size\f$*100\f$
	
	TooN::Vector<Size> delta_max; ///< Maximum distance to travel along all axes in line search

	Precision linesearch_tolerance; ///< Tolerance used to determine if the linesearch is complete. Defaults to square root of machine precision.
	Precision linesearch_epsilon;   ///< Additive term in tolerance to prevent excessive iterations if \f$x_\mathrm{optimal} = 0\f$. Known as \c ZEPS in numerical recipies. Defaults to 1e-20


	int iterations; ///< Number of iterations performed

	///Initialize the ConjugateGradient class with sensible values.
	///@param start Starting point, \e x
	///@param deriv  Function to compute \f$\nabla f(x)\f$
	///@param d Maximum allowed movement (delta) in any dimension
	template<class Deriv> ConjugateGradientOnly(const TooN::Vector<Size>& start, const Deriv& deriv, const TooN::Vector<Size>& d)
	: size(start.size()),
	  g(size),new_g(size),h(size),minus_h(size),old_g(size),old_h(size),x(start),old_x(size),delta_max(d)
	{
		init(start, deriv);
	}
	
	///Initialise given a starting point and the derivative at the starting point
	///@param start starting point
	///@param deriv derivative at starting point
	template<int S, class P, class B> void init(const TooN::Vector<Size>& start, const TooN::Vector<S, P, B>& deriv)
	{

		using std::numeric_limits;
		x = start;

		//Start with the conjugate direction aligned with
		//the gradient
		g = deriv;
		h = g;
		minus_h=-h;

		tolerance = sqrt(numeric_limits<Precision>::epsilon());
		max_iterations = size * 100;


		linesearch_tolerance =  sqrt(numeric_limits<Precision>::epsilon());
		linesearch_epsilon = 1e-20;

		iterations=0;
	}

	///Initialise given a starting point and a functor for computing derivatives
	///@param start starting point
	///@param deriv derivative computing functor
	template<class Deriv> void init(const TooN::Vector<Size>& start, const Deriv& deriv)
	{
		init(start, deriv(start));
	}

	///Perform a linesearch from the current point (x) along the current
	///conjugate vector (h).  The linesearch does not make use of values.
	///You probably do not want to call this function directly. See iterate() instead.
	///This function updates:
	/// - x
	/// - old_c
	/// - iterations
	/// Note that the conjugate direction and gradient are not updated.
	/// If bracket_minimum_forward detects a local maximum, then essentially a zero
	/// sized step is taken.
	/// @param deriv Functor returning the derivative value at a given point.
	template<class Deriv> void find_next_point(const Deriv& deriv)
	{
		iterations++;
		using std::abs;
		//Conjugate direction is -h
		//grad.-h is (should be negative)

		//We should have that f1 is negative.
		new_g = g;
		double f1 = g * minus_h;
		//If f1 is positive, then the conjugate vector points agains the
		//newly computed gradient. This can only happen if there is an error
		//in the computation of the gradient (eg we're using a stochastic method)
		//not much to do here except to stop.
		if(f1 > 0)
		{
			//Return, so try to search in a direction conjugate to the current one.
			return;
		}
		//What step size takes us up to the maximum length
		Precision max_step = HUGE_VAL;
		for(int i=0; i < minus_h.size(); i++)
			max_step = min(max_step, abs(delta_max[i]/h[i]));

		//Taking the maximum step size may result in NaNs.
		//Try the maximum step size, and seccessively reduce it.
		Precision full_max_step = max_step;
		
		for(;;)
		{
			new_g = deriv(x + max_step * minus_h);

			if(!TooN::isnan(new_g))
			{
//	cout << "new_g is NAN free :)\n";
				break;
			}
			else
				max_step /=2;

			//If the step size gets too small then 
			//return as we can't really do anything
			if(max_step < full_max_step * linesearch_tolerance)
				return;
		}

		double f2 = new_g * minus_h;


		//If f2 hasn't gone negative, then the largest allowed step is not large enough.
		//So, take a full step, then keep going in the same direction
		if(f2 < 0)
		{	
			//Take a full step
			x += max_step * minus_h;
			return;
		}

		//Now, x1 and x2 bracket a root, and find the root using bisection
		//x1 and x2 are represented by the scalars s1 and s2
		double s1 = 0;
		double s2 = max_step;
		double s_new = max_step;

		int updates[2] = {0,0};

		while(abs(s1 - s2) > abs(s1 + s2) * linesearch_tolerance + linesearch_epsilon)
		{
			if(updates[0] != updates[1] && updates[0] != 0)
			{

				//Compute the new point using false position.
				s_new = s1 + f1 * (s2 - s1)/(f1 - f2);
				new_g = deriv(x + s_new * minus_h);
				double f_new = new_g*minus_h;

				updates[0] = updates[1];

				if(f_new == 0)
					break;

				if(f_new < 0) 
				{
					s1 = s_new;
					f1 = f_new;
					updates[1] = 1;
				}
				else
				{
					s2 = s_new;
					f2 = f_new;
					updates[1] = 2;
				}
			}
			else
			{
				//Compute the new point
				
				s_new = (s1 + s2) / 2;

				new_g = deriv(x + s_new * minus_h);
				double f_new = new_g*minus_h;

				if(f_new < 0) 
				{
					s1 = s_new;
					f1 = f_new;
					updates[0] = 1;
				}
				else
				{
					s2 = s_new;
					f2 = f_new;
					updates[0] = 2;
				}
			
			}

		}

		//Update the current position and value
		//The most recent gradient computation is at s_new
		x += minus_h * s_new;

	}

	///Check to see it iteration should stop. You probably do not want to use
	///this function. See iterate() instead. This function updates nothing.
	bool finished()
	{
		using std::abs;
		return iterations > max_iterations || norm_inf(new_g) < tolerance;
	}

	///After an iteration, update the gradient and conjugate using the
	///Polak-Ribiere equations.
	///This function updates:
	///- g
	///- old_g
	///- h
	///- old_h
	///@param grad The derivatives of the function at \e x
	void update_vectors_PR(const TooN::Vector<Size>& grad)
	{
		//Update the position, gradient and conjugate directions
		old_g = g;
		old_h = h;

		g = grad;
		//Precision gamma = (g * g - oldg*g)/(oldg * oldg);
		Precision gamma = (g * g - old_g*g)/(old_g * old_g);
		h = g + gamma * old_h;
		minus_h=-h;
	}

	///Use this function to iterate over the optimization. Note that after
	///iterate returns false, g, h, old_g and old_h will not have been
	///updated.
	///This function updates:
	/// - x
	/// - old_c
	/// - iterations
	/// - g*
	/// - old_g*
	/// - h*
	/// - old_h*
	/// *'d variables not updated on the last iteration.
	///@param deriv Functor to compute derivatives at the specified point.
	///@return Whether to continue.
	template<class Deriv> bool iterate(const Deriv& deriv)
	{
		find_next_point(deriv);

		if(!finished())
		{
			update_vectors_PR(new_g);
			return 1;
		}
		else
			return 0;
	}
};


#endif

