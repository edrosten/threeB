#ifndef NUMERICAL_DERIVATIVES_H
#define NUMERICAL_DERIVATIVES_H
#include <TooN/TooN.h>

template<int S, class B, class Functor> TooN::Vector<S> debug_numerical_gradient(const Functor& f, const TooN::Vector<S,double, B>& x)
{
	using namespace TooN;
	Vector<S> grad(x.size());
	Vector<S> xh=x;
	const double h=1e-4;

	for(int i=0; i < x.size(); i++)
	{
		xh[i] += h;
		double fwd = f(xh);
		xh[i] = x[i] - h;
		double rev = f(xh);

		grad[i] = (fwd - rev) / (2*h);
		xh[i] = x[i];
	}

	return grad;
}

template<int S, class B, class Functor> TooN::Matrix<S> debug_numerical_hessian(const Functor& f, const TooN::Vector<S,double, B>& x)
{

	using namespace TooN;
	Matrix<S> hess(x.size(), x.size());
	Vector<S> xh=x;
	const double h=1e-3;

	for(int i=0; i < x.size(); i++)
		for(int j=i+1; j < x.size(); j++)
		{
			xh = x;
			xh[i] += h;
			xh[j] += h;
			double a = f(xh);

			xh = x;
			xh[i] -= h;
			xh[j] += h;
			double b = f(xh);

			xh = x;
			xh[i] += h;
			xh[j] -= h;
			double c = f(xh);

			xh = x;
			xh[i] -= h;
			xh[j] -= h;
			double d = f(xh);

			hess[i][j] = hess[j][i] = (a-b-c+d) / (4*h*h);
		}

	for(int i=0; i < x.size(); i++)
	{
		xh = x;
		double b = f(xh);

		xh[i] += h;
		double a = f(xh);
		
		xh[i] = x[i] - h;
		double c = f(xh);
		hess[i][i] = (a - 2*b + c) / (h*h);
	}

	return hess;
}


#endif
