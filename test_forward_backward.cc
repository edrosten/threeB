#include <array>
#include <tuple>
#include <TooN/TooN.h>
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

#include "forward_algorithm.h"

#undef make_tuple

using namespace std;
using namespace CVD;
using namespace tag;
using namespace TooN;

template<class C> double sumarray(const C& in)
{
	return accumulate(in.begin(), in.end(), 0., plus<double>());
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

///Observer which is the same as in hmm_test.cc
struct HmmTestObservations
{
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


int main()
{
	enableFPE();	
	cout << setprecision(10);

	HmmTestObservations Obs;
	
	double log;

	log = forward_algorithm(Obs.A, Obs.pi, Obs, Obs.O);

	cout << "From forward_algorithm: " << log << endl;

	vector<array<double, 2> > delta, epsilon;

	tie(delta, epsilon) = forward_backward_algorithm(Obs.A, Obs.pi, Obs, Obs.O);

	for(unsigned int i=0; i < delta.size(); i++)
	{
		double sum=0;
		for(unsigned int j=0; j < delta[j].size(); j++)
			sum += exp(delta[i][j] + epsilon[i][j]);
		cout << ln(sum) << endl;
	}
	array<double, 2> zero;
	zero.assign(0);

	vector<array<double, 2> > samples(delta.size(), zero);

	RngDrand48 rng;

	for(int i=0; i < 100000; i++)
	{
		vector<int> s = backward_sampling<2, int>(Obs.A, delta, rng, 1e10);

		for(unsigned int j=0; j < s.size(); j++)
			samples[j][s[j]]++;
	}

	cout << setprecision(5);
	for(unsigned int i=0; i < delta.size(); i++)
	{
		double hi=0;
		for(unsigned int j=0; j < delta[j].size(); j++)
			hi = max(delta[i][j] + epsilon[i][j], hi);

		double sum=0;
		for(unsigned int j=0; j < delta[j].size(); j++)
			sum += exp(delta[i][j] + epsilon[i][j] - hi);

		double lnp = ln(sum) + hi;
		
		for(unsigned int j=0; j < delta[j].size(); j++)
			cout << setw(10) << exp(delta[i][j] + epsilon[i][j] - lnp) << " ";

		cout << "    ";
		double s = sumarray(samples[i]);
		
		for(unsigned int j=0; j < delta[j].size(); j++)
			cout << setw(10) << samples[i][j]/s << " ";

		cout << " " << lnp << endl;
	}

}
