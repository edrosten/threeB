#ifndef MT19937_H
#define MT19937_H

#include "randomc.h"
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>

///Useful wrapper for MT19937 random number generator class.
///@ingroup gUtility
struct MT19937
{
	///Null struct thrown if attempting to load
	///state from stream yields a parse error.
	struct ParseError{};
	
	///Underlying RNG.
	CRandomMersenne rng;

	public:
		MT19937()
		:rng(0)
		{}
		
		///Seed state with a simple RNG
		///@param seed
		void simple_seed(int seed)
		{
			rng.RandomInit(seed);
		}

		///Duplicate RNG state
		///@param r RNG to duplicate
		void copy_state(const MT19937& r)
		{
			rng = r.rng;
		}
		
		///Generate a double
		double operator()()
		{
			return rng.Random();
		}

		///Generate an int
		uint32_t rand_int()
		{
			return rng.BRandom();
		}
		
		///Generate a Gaussian variate
		double gaussian()
		{
			double x1, x2, w;
			do {
				x1 = 2.0 * (*this)() - 1.0;
				x2 = 2.0 * (*this)() - 1.0;
				w = x1 * x1 + x2 * x2;
			} while ( w >= 1.0 );

			w = std::sqrt( (-2.0 * std::log( w ) ) / w );
			return x1 * w;
			//spare so we don't have to save that one extra bit of state y2 = x2 * w;
		}
		
		///Serialise state
		///@param o Stream to serialise to
		void write(std::ostream& o)
		{
			using namespace std;
			char f = o.fill();
			ios_base::fmtflags fl = o.flags();
			o << "MT19937 " << hex << setfill('0') << setw(3) << rng.get_index();	
			for(int i=0; i < MERS_N; i++)
				o << " " << hex << setw(8) << rng.get_state()[i];
			
			o << setfill(f) << setiosflags(fl);
		}

		///De serialise state
		///param is Stream to de-serialise from
		void read(std::istream& is)
		{
			using namespace std;

			string ls;
			getline(is, ls);
			if(ls.size() != 5627)
			{
				cerr << "MT19937: Expected string of length 5627. Got " << ls.size() << endl;
				throw ParseError();
			}

			istringstream l(ls);
			
			string s;
			uint32_t i;

			l >> s;
			
			if(s != "MT19937")
			{	
				cerr << "MT19937: Expected MT19937. Got " << s << endl;
				throw ParseError();
			}		

			for(int n=0; n < MERS_N + 1; n++)
			{
				l >> hex >> i;
				if(l.bad())
				{
					cerr << "MT19937: Expected hex number. Got ";
					if(l.eof())
						cerr << "EOF" << endl;
					else
					{
						cerr << l.get() << endl;
					}

					throw ParseError();
				}

				if(n==0)
					rng.get_index() = i;
				else
					rng.get_state()[n-1]=i;

			}

		}
	private:
		/// Disallow copying, since one almost never wants to do this.
		/// Copying has to be explicit via MT19937::copy_state().
		MT19937(const MT19937&);

};


#endif
