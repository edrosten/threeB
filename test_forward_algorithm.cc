#include <array>
#include <tuple>
#include <TooN/TooN.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <tag/stdpp.h>
#include <cvd/cpu_hacks.h>

#include "forward_algorithm.h"
#include "numerical_derivatives.h"

#undef make_tuple

using namespace std;
using namespace CVD;
using namespace tag;
using namespace TooN;

///Observer which is the same as in hmm_test.cc
struct HmmTestObservations
{
	static const int NumParameters = 6;
	//Pobability of emmiting symbols in a given state.
	Vector<6> parameters, real_parameters;
	Matrix<2, 3, double, Reference::RowMajor> B;
	Vector<2> pi;
	Matrix<2> A;
	vector<int> O, Q;

	HmmTestObservations()
	:real_parameters(makeVector(.6,.2,.2,.2,.1,.7)),B(&parameters[0])
	{
		//B[0] = makeVector(.6,.2,.2);  //Probabilities of symbols being emitted in state 0
		//B[1] = makeVector(.2,.1,.7);  //Probabilities of symbols being emitted in state 1

		parameters=real_parameters;

		srand48(0);
		//Transition probabilities.

		//Row  = from col = to
		A[0] = makeVector(.9, .1);
		A[1] = makeVector(.2, .8);
		
		//Initial state probabilities
		//Note numerical derivatives for an intermadiate result
		//will appear incorrect if any of these are exactly 0
		pi = makeVector(.9, .1);

		Q = run_hmm(A, pi, 100);
		O = make_observations(B, Q);

	}

	double log(int state, int observation) const
	{
		assert(state == 0 || state == 1);
		assert(observation >=0 && observation < 3);
		return ::log(B[state][observation]);
	}

	Vector<NumParameters> diff_log(int state, int observation) const
	{
		assert(state == 0 || state == 1);
		assert(observation >=0 && observation < 3);
		Vector<6>  ret = Zeros;

		ret[observation + state*3] = 1/B[state][observation];

		return ret;
	}

	Matrix<NumParameters> hess_log(int state, int observation) const
	{
		assert(state == 0 || state == 1);
		assert(observation >=0 && observation < 3);
		return -diff_log(state, observation).as_col() * diff_log(state, observation).as_row();
	}

	tuple<double, Vector<NumParameters>, Matrix<NumParameters> > log_diff_hess(int state, int observation) const
	{
		return make_tuple(log(state, observation), diff_log(state, observation), hess_log(state, observation));
	}

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

	template<int States, int Outputs, class Base> vector<int> make_observations(const Matrix<States, Outputs, double, Base>& B, const vector<int>& Q)
	{
		vector<int> O;
		for(unsigned int i=0; i < Q.size(); i++)
			O.push_back(select_random_element(B[Q[i]]));

		return O;
	}
};

struct EvalHmmTestObservations
{
	double operator()(const Vector<6>& x) const
	{
		HmmTestObservations h;
		h.parameters = x;
		return forward_algorithm(h.A, h.pi, h, h.O);
	}

};



int main()
{
	enableFPE();	
	cout << setprecision(10);

	HmmTestObservations Obs;
	
	double log;
	Vector<6> deriv;
	Matrix<6> hess;

	tie(log, deriv, hess) = forward_algorithm_hessian(Obs.A, Obs.pi, Obs, Obs.O);

	cout << "Correct answer = -101.4079256\n";
	cout << "------------- Computed values ---------------\n";
	cout << log << endl << endl;
	cout << deriv << endl << endl;
	cout << hess << endl;

	cout << "------------- Numerical Derivatives---------------\n";
	cout << debug_numerical_gradient(EvalHmmTestObservations(), makeVector(.6,.2,.2,.2,.1,.7)) << endl << endl;
	cout << debug_numerical_hessian(EvalHmmTestObservations(), makeVector(.6, .2, .2, .2, .1, .7))  << endl;
}
