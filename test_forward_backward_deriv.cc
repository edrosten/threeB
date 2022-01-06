#include <array>
#include <tuple>
#include <TooN/TooN.h>
#include <TooN/functions/derivatives.h>
#include <algorithm>
#include <numeric>
#include <functional>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <tag/stdpp.h>
#include <cvd/cpu_hacks.h>
#include <cvd/random.h>

#include "forward_algorithm.h"
#include "utility.h"

#undef make_tuple

using namespace std;
using namespace CVD;
using namespace tag;
using namespace TooN;


template<int I, class Base> int select_random_element(const Vector<I, double, Base>& v)
{
	double total=0, choice = drand48();

	for(int i=0; i < v.size(); i++)
	{
		total += v[i];	
		if(choice <= total)
			return i;
	}
	return v.size()-1;
}


template<int States> vector<int> run_hmm(Matrix<States> A, Vector<States> pi, int n)
{
	int state = select_random_element(pi);
	vector<int> states;
	for(int i=0 ;i<n; i++)
	{
		states.push_back(state);
		state = select_random_element(A[state]);
	}
	return states;
}

template<int States> vector<double> generate_gaussian_observations(const vector<int>& states, const Vector<States*2>& params)
{
	vector<double> obs;

	for(unsigned int i=0; i < states.size(); i++)
		obs.push_back(rand_g() * params[states[i]*2+1] + params[states[i]*2]);
	return obs;
}

template<int States> double log_prob(double obs, int state, const Vector<States*2>& params)
{
	double mu = params[state*2];
	double sigma = params[state*2+1];
	return -sq(obs - mu) / (2*sq(sigma)) - .5*log(2*M_PI) - log(sigma);
}

template<int States> Vector<States*2> prob_diff(double obs, int state, const Vector<States*2>& params)
{
	double mu = params[state*2];
	double sigma = params[state*2+1];
		
	Vector<States*2> ret = Zeros;
	
	//d foo / d mu
	ret[state*2+0] = (obs - mu) / sq(sigma); 
	ret[state*2+1] = sq(obs - mu) / (sigma*sigma*sigma) - 1/sigma;

	return ret;
}

template<int States> struct BjOt
{
	static const int NumParameters = States*2;
	Vector<States*2> p;

	BjOt(const Vector<States*2>& v)
	:p(v)
	{}

	double log(int state, double obs) const
	{
		return log_prob<States>(obs, state, p);
	}

	Vector<States*2> diff_log(int state, double obs) const
	{
		return prob_diff<States>(obs, state, p);
	}
	
	Matrix<States*2> hess_log(int, double) const
	{
		return Zeros;
	}
};


struct EvalHmm{
	vector<double> obs;
	Matrix<2> A;
	Vector<2> pi;

	double operator()(const Vector<4>& p) const
	{
		return forward_algorithm(A, pi, BjOt<2>(p), obs);
	}
};


int main()
{
	enableFPE();	
	cout << setprecision(10);

	Matrix<2> A = Data(.5, .5, .5, .5);
	Vector<2> pi = makeVector(.5,.5);
	Vector<4> params = makeVector(0, 10, 0, 1);

	vector<int>  states = run_hmm(A, pi, 10);
	vector<double> O = generate_gaussian_observations<2>(states, params);


	BjOt<2> B(params);


	cout << forward_algorithm(A, pi, B, O) << endl;
	cout << "Exact with forward algo = " << forward_algorithm_deriv(A, pi, B, O).second << endl;

	EvalHmm e;
	e.A = A;
	e.pi = pi;
	e.obs = O;
	cout << "Approx with finite diff = " << numerical_gradient(e, params) << endl;

	Vector<4> sum_diff = Zeros;
	vector<array<double, 2> > delta, epsilon;
	tie(delta, epsilon) = forward_backward_algorithm(A, pi, B, O);
	
	int N = 100000;
	for(int i=0; i < N; i++)
	{
		vector<int> s = backward_sampling<2>(A, delta);

		Vector<4> E = Zeros;
		for(unsigned int j=0; j < s.size(); j++)
		{
			E += prob_diff<2>(O[j], s[j], params);
		}
		//cout << E << endl;
		sum_diff+=E;
	}

	cout << "Approx with ffbs sampling " << sum_diff / N << endl;


}

