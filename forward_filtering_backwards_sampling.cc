#include <array>
#include <TooN/TooN.h>
#include <vector>
#include <cmath>

#include "forward_algorithm.h"


template<int States, class Btype, class Otype> void forward_backward_sample(TooN::Matrix<States> A, TooN::Vector<States> pi, const Btype& B, const std::vector<Otype>& O)
{
	using namespace TooN;
	using namespace std;

	static const int M = Btype::NumParameters;
	int states = pi.size();
	
	//delta[j][i] = delta_t(i)
	vector<array<double, States> > delta(O.size());

	//Initialization: Eqn 19, P 262
	//Set initial partial log probabilities:
	for(int i=0; i < states; i++)
		delta[0][i] = ln(pi[i]) + B.log(i, O[0]);

	//Perform the recursion: Eqn 20, P262
	//Note, use T and T-1. Rather than  T+1 and T.
	for(unsigned int t=1; t < O.size(); t++)
	{	
		for(int j=0; j < states; j++)
		{
			double Ztj = -HUGE_VAL; //This is Z_t(j)
			for(int i=0; i < states; i++)
				Ztj = max(Ztj, delta[t-1][i] + ln(A[i][j]));

			double sum=0;
			for(int i=0; i < states; i++)
				sum += exp(delta[t-1][i] + ln(A[i][j]) - Ztj);

			delta[t][j] = B.log(j, O[t]) + Ztj + ln(sum);
		}
	}

	//Compute the log prob using normalization
	double Z = -HUGE_VAL;
	for(int i=0; i < states; i++)
		Z = max(Z, delta.back()[i]);

	double sum =0;
	for(int i=0; i < states; i++)
		sum += exp(delta.back()[i] - Z);

	double log_prob = Z  + ln(sum);

	//Compute the differential of the log
	Vector<M> diff_log = Zeros;
	//Compute the differential of the log using normalization
	//The convenient normalizer is ln P(O|lambda) which makes the bottom 1.
	for(int i=0; compute_deriv && i < states; i++)
		diff_log += exp(delta.back()[i] - log_prob)*diff_delta.back()[i];

	Matrix<M> hess_log = Zeros;
	//Compute the hessian of the log using normalization
	//The convenient normalizer is ln P(O|lambda) which makes the bottom 1.
	for(int i=0; compute_hessian && i < states; i++)
		hess_log += exp(delta.back()[i] - log_prob) * (hess_delta.back()[i] + diff_delta.back()[i].as_col() * diff_delta.back()[i].as_row());

	hess_log -= diff_log.as_col() * diff_log.as_row();

	//Compute the differential of the Hessian
	return  make_tuple(log_prob, diff_log, hess_log);
}


