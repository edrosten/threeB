#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "mt19937.h"
#include "poisson.h"

#include <vector>

using namespace std;

namespace {
void set_mt19937(void* r, unsigned long int seed)
{
	MT19937& rng = *((MT19937*)r);
	rng.simple_seed(seed);
}

unsigned long int sample_int_mt19937(void* r)
{
	MT19937& rng = *((MT19937*)r);
	return rng.rand_int();
}

double sample_double_mt19937(void* r)
{
	MT19937& rng = *((MT19937*)r);
	return rng();
}

static const gsl_rng_type rng_local_19937 = 
{
	"local MT19937",
	0xffffffff,
	0,
	0,
	set_mt19937,
	sample_int_mt19937,
	sample_double_mt19937
};


static double fact(int n)
{
	double f=1;
	for(int i=2; i <= n; i++)
		f*=i;
	return f;
}
}

unsigned int poisson(double mu, MT19937& rng)
{
	gsl_rng grng = { &rng_local_19937, &rng};

	return gsl_ran_poisson (&grng, mu);
}

/*int main()
{
	MT19937 rng;
	vector<int> hist(200);
	int sum=0;
	double mu=100.5;

	for(int i=0; i < 10000000; i++)
	{
		unsigned int n = fish(mu, rng);

		if(n < hist.size())
		{
			hist[n]++;
			sum++;
		}
	}

	for(unsigned int i=0; i < hist.size(); i++)
	{
		cout << i << " " << hist[i]*1.0/sum << " " << (pow(mu, i)/fact(i) * exp(-mu)) << endl;
	}

}*/
