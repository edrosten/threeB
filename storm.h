#ifndef INC_STORM_H
#define INC_STORM_H

#include <TooN/TooN.h>
#include <cvd/image.h>
#include <utility>
#include <tuple>
#include "utility.h"


/**See spot_shape()
@param x \f$\Vec{x}\f$
@param phi \f$\Vec{\phi}\f$
@return \f$s(\Vec{x}, \Vec{\phi}) \f$
@ingroup gStorm
*/
template<class B> double spot_shape_s(const TooN::Vector<2>& x, const TooN::Vector<4, double, B>& phi)
{
	return -norm_sq(x - phi.template slice<2,2>()) / (2*phi[1]*phi[1]);
}

/** Compute the spot shape and its derivative with respect to posision. See also spot_shape()
@param x \f$\Vec{x}\f$
@param phi \f$\Vec{\phi}\f$
@ingroup gStorm
*/
template<class B> std::pair<double, TooN::Vector<4> > spot_shape_diff_position(const TooN::Vector<2>& x, const TooN::Vector<4, double, B>& phi)
{
	using namespace TooN;

	double s = spot_shape_s(x, phi);
	double r_2_pi = sqrt(2*M_PI);

	double prob = exp(s) * phi[0]/(phi[1]*r_2_pi);
	
	Vector<4> deriv = (exp(s) / (phi[1]*r_2_pi)) * 
	                       makeVector(1, 
						              -phi[0] * (1 + 2*s)/phi[1], 
									  (x[0] - phi[2])*(phi[0]/sq(phi[1])), 
									  (x[1] - phi[3])*(phi[0]/sq(phi[1])));
	return std::make_pair(prob, deriv);
}

/** Compute the spot shape and its Hessian with respect to posision. See also spot_shape()
@param x \f$\Vec{x}\f$
@param phi \f$\Vec{\phi}\f$
@ingroup gStorm
*/
template<class B> std::tuple<double, TooN::Vector<4>, TooN::Matrix<4> > spot_shape_hess_position(const TooN::Vector<2>& x, const TooN::Vector<4, double, B>& phi)
{
	using namespace TooN;
	using namespace std;

	double s = spot_shape_s(x, phi);
	double r_2_pi = sqrt(2*M_PI);

	double es = exp(s);

	double prob = es * phi[0]/(phi[1]*r_2_pi);
	
	Vector<4> deriv = (es / (phi[1]*r_2_pi)) * 
	                       makeVector(1, 
						              -phi[0] * (1 + 2*s)/phi[1], 
									  (x[0] - phi[2])*(phi[0]/sq(phi[1])), 
									  (x[1] - phi[3])*(phi[0]/sq(phi[1])));

	Matrix<4> hess;
	hess[0][0] = 0;

	hess[0][1] = -es*(1+2*s) / (phi[1] * phi[1] * r_2_pi);
	hess[1][0] = hess[0][1];

	hess[0][2] = es * (x[0] - phi[2]) / (pow(phi[1], 3)*r_2_pi);
	hess[2][0] = es * (x[0] - phi[2]) / (pow(phi[1], 3)*r_2_pi);

	hess[0][3] = es * (x[1] - phi[3]) / (pow(phi[1], 3)*r_2_pi);
	hess[3][0] = es * (x[1] - phi[3]) / (pow(phi[1], 3)*r_2_pi);

	hess[1][1] = 2*phi[0]*es*(1 + 5*s + 2*s*s) / ( pow(phi[1], 3) * r_2_pi);

	hess[1][2] = -phi[0] * es * (3 + 2*s) * (x[0] - phi[2]) / (pow(phi[1], 4) * r_2_pi);
	hess[1][3] = -phi[0] * es * (3 + 2*s) * (x[1] - phi[3]) / (pow(phi[1], 4) * r_2_pi);

	hess[2][1] = hess[1][2];
	hess[3][1] = hess[1][3];

	hess[2][2] = phi[0] * es * (sq(x[0] - phi[2]) - sq(phi[1])) / (r_2_pi * pow(phi[1], 5));
	hess[3][3] = phi[0] * es * (sq(x[1] - phi[3]) - sq(phi[1])) / (r_2_pi * pow(phi[1], 5));
	
	hess[2][3] = phi[0] * es * (x[0] - phi[2])*(x[1] - phi[3]) / (r_2_pi * pow(phi[1], 5));
	hess[3][2] = hess[2][3];


	return make_tuple(prob, deriv, hess);
}

/**Value of the spot, given the parameters and input location.
The spot is described by the following formula:
\f[
	\mu(\Vec{x}, \Vec{\phi}) = \frac{\phi_1}{\phi_2\sqrt(2\pi)} e^s,
\f]
where
\f[
	s = -\frac{(x_1 - \phi_3)^2 + (x_2 - \phi_4)^2}{2\phi_2^2}.
\f]
This describes a generic blobby spot function of a variable size. The light output
can be tuned by varying \f$\phi_1\f$, and the level of blur can be changed independently
by varying \f$\phi_2\f$. The derivative is:
\f{eqnarray}{
	\frac{\partial \mu}{\partial \phi_1} &=& \frac{1}{\phi_2\sqrt{2\pi}}e^s\\
	\frac{\partial \mu}{\partial \phi_2} &=& -\frac{\phi_1}{\phi_2^2\sqrt{2\pi}}e^s(1 + 2s)
\f}
And the hessian is:
\f{eqnarray}{
	\frac{\partial^2 \mu}{\partial \phi_1^2} &=& 0\\
	\frac{\partial^2 \mu}{\partial\phi_1 \partial \phi_2} &=& -\frac{1}{\phi_2^2\sqrt{2\pi}}e^s(1 + 2s)\\
	\frac{\partial^2 \mu}{\partial \phi_2^2} &=& \frac{2\phi_1}{\phi_2^3\sqrt{2\pi}}e^s(1 + 5s + 2s^2)
\f}
@param x \f$\Vec{x}\f$
@param phi \f$\Vec{\phi}\f$
@return \f$\mu(\Vec{x}, \Vec{\phi}) \f$
@ingroup gStorm
*/
template<class B> std::tuple<double, TooN::Vector<2>, TooN::Matrix<2> > spot_shape_hess(const TooN::Vector<2>& x, const TooN::Vector<4, double, B>& phi)
{
	double s = spot_shape_s(x, phi);
	double r_2_pi = sqrt(2*M_PI);

	double prob = exp(s) * phi[0]/(phi[1]*r_2_pi);
	TooN::Vector<2> deriv = (exp(s) / (phi[1]*r_2_pi)) * TooN::makeVector(1, -phi[0] * (1 + 2*s)/phi[1]);
	TooN::Matrix<2> hess;

	hess[0][0] = 0;
	hess[0][1] = -exp(s)*(1+2*s) / (phi[1] * phi[1] * r_2_pi);
	hess[1][0] = hess[0][1];
	hess[1][1] = 2*phi[0]*exp(s)*(1 + 5*s + 2*s*s) / ( pow(phi[1], 3) * r_2_pi);

	return std::make_tuple(prob, deriv, hess);
}
/** see spot_shape_hess()
@param x \f$\Vec{x}\f$
@param phi \f$\Vec{\phi}\f$
@return \f$\mu(\Vec{x}, \Vec{\phi}) \f$
@ingroup gStorm
*/
template<class B> std::pair<double, TooN::Vector<2> > spot_shape_diff(const TooN::Vector<2>& x, const TooN::Vector<4, double, B>& phi)
{
	double s = spot_shape_s(x, phi);
	double r_2_pi = sqrt(2*M_PI);

	double prob = exp(s) * phi[0]/(phi[1]*r_2_pi);
	TooN::Vector<2> deriv = (exp(s) / (phi[1]*r_2_pi)) * TooN::makeVector(1, -phi[0] * (1 + 2*s)/phi[1]);
	return std::make_pair(prob, deriv);
}

/** see spot_shape_hess()
@param x \f$\Vec{x}\f$
@param phi \f$\Vec{\phi}\f$
@return \f$\mu(\Vec{x}, \Vec{\phi}) \f$
@ingroup gStorm
*/
template<class B> double spot_shape(const TooN::Vector<2>& x, const TooN::Vector<4, double, B>& phi)
{
	double s = spot_shape_s(x, phi);
	double r_2_pi = sqrt(2*M_PI);
	
	// FIXME FIXME FIXME and don't forget to fix the HESSIAN AND DERIVATIVE
	// Should be:              1/(2 pi s^2)  for two dimensions
	//                             vvvvvvvvvvvvv    http://lol.i.trollyou.com/
	double prob = exp(s) * phi[0]/(phi[1]*r_2_pi);


	return prob;
}

/**Find the log probability of an image patch, 
assuming zero mean and the given variance, and no spot present.
See also log_probability_spot()
@param im Image
@param variance variance
@returns The log probability
*/
inline double log_probability_no_spot(const CVD::SubImage<float>& im, double variance)
{
	double logprob_part=0;
	for(int y=0; y < im.size().y; y++)
		for(int x=0; x < im.size().x; x++)
			logprob_part -= im[y][x] * im[y][x];
	return logprob_part/(2*variance) - im.size().area() * log(2*M_PI*variance)/2;

}

/**Find the log probability of an image patch, assuming zero base-line mean and
the given variance. This function makes use of the spot shape. It is assumed that the
centre pixel of the image is at 0,0. Since the noise is Gaussian:
\f{eqnarray}{
	P(\text{image})    &=& \prod_{\Vec{x} \in \text{pixels}} \frac{1}{\sqrt{2\pi\sigma^2}}e^{-\frac{(I(\Vec{x}) - \mu(\Vec{x}, \Vec{\phi}))^2}{2\sigma^2}} \\
	\ln P(\text{image}) &=& \sum_{\Vec{x} \in \text{pixels}}  -\frac{(I(\Vec{x}) - \mu(\Vec{x}, \Vec{\phi}))^2}{2\sigma^2} - \frac{N}{2} \ln {2 \pi \sigma^2},
\f}
where \e I is the image, and \e N is the number of pixels. See also ::log_probability_no_spot and \f$\mu\f$ (::spot_shape).
The derivatives are:
\f{eqnarray}{
	\frac{\partial \ln P(I)}{\partial \phi_0} &=& \frac{1}{\sigma^2} \sum_{\Vec{x}}(I_{\Vec{x}} - \mu(\Vec{x},\Vec{\phi}))
                                                           \frac{\partial}{\partial \phi_0}\mu(\Vec{x}, \Vec{\phi})\\
	\frac{\partial^2 \ln P(I)}{\partial \phi_0 \partial \phi_1} &=&
	          \frac{1}{\sigma^2} \sum_{\Vec{x}}(I_{\Vec{x}} - \mu(\Vec{x},\Vec{\phi}))
                               \frac{\partial^2}{\partial \phi_0 \partial \phi_1}\mu(\Vec{x}, \Vec{\phi}) - 
							   \frac{\partial}{\partial \phi_0}\mu(\Vec{x},\Vec{\phi})
							   \frac{\partial}{\partial \phi_1}\mu(\Vec{x},\Vec{\phi})
\f}
@ingroup gStorm
@param im Image
@param variance \f$\sigma^2\f$
@param spot_parameters \f$\Vec{\phi}\f$
@returns The log probability
*/
template<class Base> std::tuple<double, TooN::Vector<2>, TooN::Matrix<2> > log_probability_spot_hess(const CVD::SubImage<float>& im, double variance, const TooN::Vector<4, double, Base>& spot_parameters)
{
	using namespace TooN;
	using namespace std;

	//-1 because if the image is 3x3, ie 0,1,2 then 1,1 is the centre.
	//If it is 2x2, ie 0,1 then .5,.5 is the centre
	Vector<2> centre = makeVector((im.size().x-1) / 2.0, (im.size().y-1) / 2.0);

	double logprob_part=0;
	Vector<2> diff = Zeros;
	Matrix<2> hess = Zeros;
	for(int y=0; y < im.size().y; y++)
		for(int x=0; x < im.size().x; x++)
		{
			Vector<2> d = TooN::makeVector(x, y) - centre;

			double mu;
			Vector<2> diff_mu;
			Matrix<2> hess_mu;
			tie(mu, diff_mu, hess_mu) = spot_shape_hess(d, spot_parameters);

			double e = im[y][x] - mu;

			logprob_part += -sq(e);
			diff         += diff_mu * e;
			hess         += e * hess_mu - diff_mu.as_col() * diff_mu.as_row();
		}
	return make_tuple(   logprob_part / (2*variance) - im.size().area() * log(2*M_PI*variance)/2,
						diff / variance,
						hess / variance);
}

/** See log_probability_spot_hess
@ingroup gStorm
@param im Image
@param variance \f$\sigma^2\f$
@param spot_parameters \f$\Vec{\phi}\f$
@returns The log probability
*/
template<class Base> std::pair<double, TooN::Vector<2> > log_probability_spot_diff(const CVD::SubImage<float>& im, double variance, const TooN::Vector<4, double, Base>& spot_parameters)
{
	using namespace TooN;
	using namespace std;
	//-1 because if the image is 3x3, ie 0,1,2 then 1,1 is the centre.
	//If it is 2x2, ie 0,1 then .5,.5 is the centre
	Vector<2> centre = makeVector((im.size().x-1) / 2.0, (im.size().y-1) / 2.0);

	double logprob_part=0;
	Vector<2> diff = Zeros;
	for(int y=0; y < im.size().y; y++)
		for(int x=0; x < im.size().x; x++)
		{
			Vector<2> d = makeVector(x, y) - centre;

			double mu;
			Vector<2> diff_mu;
			tie(mu, diff_mu) = spot_shape_diff(d, spot_parameters);

			double e = im[y][x] - mu;

			logprob_part += -sq(e);
			diff         += diff_mu * e;
		}
	return make_pair(logprob_part / (2*variance) - im.size().area() * log(2*M_PI*variance)/2, diff / variance);
}

/** See log_probability_spot_hess
@ingroup gStorm
@param im Image
@param variance \f$\sigma^2\f$
@param spot_parameters \f$\Vec{\phi}\f$
@returns The log probability
*/
template<class Base> double log_probability_spot(const CVD::SubImage<float>& im, double variance, const TooN::Vector<4, double, Base>& spot_parameters)
{
	//-1 because if the image is 3x3, ie 0,1,2 then 1,1 is the centre.
	//If it is 2x2, ie 0,1 then .5,.5 is the centre
	TooN::Vector<2> centre = TooN::makeVector((im.size().x-1) / 2.0, (im.size().y-1) / 2.0);

	double logprob_part=0;
	for(int y=0; y < im.size().y; y++)
		for(int x=0; x < im.size().x; x++)
		{
			TooN::Vector<2> d = TooN::makeVector(x, y) - centre;

			double mu = spot_shape(d, spot_parameters);

			double e = im[y][x] - mu;

			logprob_part += -sq(e);
		}
	return logprob_part / (2*variance) - im.size().area() * log(2*M_PI*variance)/2;
}

/**Compute the standard deviation of a log-normal distribution.

See log_normal().
\f{equation}
	\mathrm{Var}[P(x)] = (e^(\sigma^2)-1)e^(2*\mu+\sigma^2)
\f}
@param sigma \f$ \sigma\f$
@param mu \f$ \mu\f$
@returns The standard deviation
@ingroup gStorm
*/
inline double log_normal_std(double mu, double sigma)
{
	return sqrt((exp(sq(sigma)) - 1) * exp(2*mu + sq(sigma)));
}

/**Compute the mode of a log-normal distribution.

See log_normal().
\f{equation}
	\mathrm{Mode}[P(x)] = e^(\mu-\sigma^2)
\f}
@param sigma \f$ \sigma\f$
@param mu \f$ \mu\f$
@returns The mode
@ingroup gStorm
*/
inline double log_normal_mode(double mu, double sigma)
{
	return exp(mu - sigma * sigma);
}

/**Log-normal distribution. This is given by:
\f{eqnarray}{
	P(x)     &=& \frac{1}{x\sigma\sqrt{2\pi}} e^{-\frac{(\ln x - \mu)^2}{s\sigma^2}}\\
    \ln P(x) &=& -\frac{(\ln x - \mu)^2}{s\sigma^2} - \ln x - \ln\sigma\sqrt{2\pi}.
\f}
@param x \e x
@param mu \f$\mu\f$
@param sigma \f$\sigma\f$
@ingroup gStorm
*/
inline double log_log_normal(double x, double mu, double sigma)
{
	return -sq(ln(x) - mu) / (2*sq(sigma)) - ln(x) - ln(sigma * sqrt(2*M_PI));
}

/**Derivative of the log of the log-normal distribution:
\f[
	\frac{\partial \ln P(x)}{\partial x} = -\frac{1}{x}\left(1 + \frac{\ln x - \mu}{\sigma^2}\right).
\f]
@param x \e x
@param mu \f$\mu\f$
@param sigma \f$\sigma\f$
@ingroup gStorm
*/
inline double diff_log_log_normal(double x, double mu, double sigma)
{
	return -(1 + (ln(x) - mu)/sq(sigma)) / x;
}


/**Second derivative of the log of the log-normal distribution:
\f[
	\frac{\partial^2 \ln P(x)}{\partial x^2} = \frac{1}{x^2}\left(1 + \frac{\ln x - \mu}{\sigma^2} - \frac{1}{\sigma^2}\right).
\f]
@param x \e x
@param mu \f$\mu\f$
@param sigma \f$\sigma\f$
@ingroup gStorm
*/
inline double hess_log_log_normal(double x, double mu, double sigma)
{
	return (1 + (ln(x) - mu - 1)/sq(sigma)) / sq(x);
}


#endif
