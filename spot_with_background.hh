
///This class compute the log-diff-hess probability of a spot, given an
///image patch and background due to existing spots.
struct SWBG_NAME
{
	static const int NumParameters=4;

	vector<pair<double, double> > log_prob;
	vector<Vector<4> > diff_log_prob;
	vector<Matrix<4> > hess_log_prob;

	double    get_val (double d){return d;}
	Vector<4> get_diff(double){return Zeros;}
	Matrix<4> get_hess(double){return Zeros;}

	double    get_val (const pair<double, Vector<4> >& d){return d.first;}
	Vector<4> get_diff(const pair<double, Vector<4> >& d){return d.second;}
	Matrix<4> get_hess(const pair<double, Vector<4> >&){return Zeros;}

	double    get_val (const tuple<double, Vector<4>, Matrix<4> >& d){return get<0>(d);}
	Vector<4> get_diff(const tuple<double, Vector<4>, Matrix<4> >& d){return get<1>(d);}
	Matrix<4> get_hess(const tuple<double, Vector<4>, Matrix<4> >& d){return get<2>(d);}

	template<class C> bool type_has_hess(const C&)                                   { return 0;}
	                  bool type_has_hess(const tuple<double, Vector<4>, Matrix<4> >&){ return 1;}

	template<class C> bool type_has_diff(const C&)                                   { return 0;}
	                  bool type_has_diff(const pair<double, Vector<4> > &)           { return 1;}
	                  bool type_has_diff(const tuple<double, Vector<4>, Matrix<4> >&){ return 1;}
	
#ifdef SWBG_HAVE_DRIFT
	#define SWBG_SPOT_INTENSITIES vector<vector<Input> >
#else
	#define SWBG_SPOT_INTENSITIES vector<Input>
#endif

	template<class Input>
#ifdef SWBG_HAVE_MASK
	SWBG_NAME(const vector<vector<double> >& sample_intensities, const SWBG_SPOT_INTENSITIES& spot_intensities, const vector<vector<double> >& pixel_intensities, const double variance, const vector<int>& mask)
#else
	SWBG_NAME(const vector<vector<double> >& sample_intensities, const SWBG_SPOT_INTENSITIES& spot_intensities, const vector<vector<double> >& pixel_intensities, const double variance)
#endif

	{
		//The constructor computes log_probability_spot and log_probability_no_spot
		
		//Indexes are:
		//  sample_intensities[frame][pixel]
		//  spot_intensities[[frame]][pixel]
		//  pixel_intensities[frame][pixel]


#ifdef SWBG_HAVE_MASK
		const unsigned int n_pix = mask.size();
		#define MASK(X) mask[X]
#else
		const unsigned int n_pix = sample_intensities[0].size();
		#define MASK(X) X
#endif
		const int n_frames = sample_intensities.size();

		assert(sample_intensities.size() == pixel_intensities.size());
		assert_same_size(sample_intensities);
		assert_same_size(pixel_intensities);
#ifdef SWBG_HAVE_DRIFT
		assert_same_size(spot_intensities);
#endif
		
		log_prob.resize(sample_intensities.size());

#ifdef SWBG_HAVE_DRIFT
		const bool have_diff = type_has_diff(spot_intensities[0][0]);
		const bool have_hess = type_has_hess(spot_intensities[0][0]);
#else
		const bool have_diff = type_has_diff(spot_intensities[0]);
		const bool have_hess = type_has_hess(spot_intensities[0]);
#endif

		if(have_diff)
			diff_log_prob.resize(log_prob.size());

		if(have_hess)
			hess_log_prob.resize(log_prob.size());
		
		//Compute all probabilities
		for(int frame=0; frame < n_frames; frame++)
		{
			//assert(pixel_intensities[frame].size() == n_pix);
			//assert(sample_intensities[frame].size() == n_pix);

			double log_prob_off=0, log_prob_on = 0;
			Vector<4> diff = Zeros;
			Matrix<4> hess = Zeros;


			for(unsigned int p=0; p < n_pix; p++)
			{

#ifdef SWBG_HAVE_DRIFT
				const int n_steps = spot_intensities.size();
				int s = frame_to_step(n_frames, n_steps, frame);
	#define spot_intensities spot_intensities[s] 
#else
	#define spot_intensities spot_intensities
#endif

				double e = pixel_intensities[frame][MASK(p)] - (sample_intensities[frame][MASK(p)] + get_val(spot_intensities[MASK(p)]));


				log_prob_off -= sq(pixel_intensities[frame][MASK(p)] - sample_intensities[frame][MASK(p)]);
				log_prob_on  -= sq(e);
				
				if(have_diff)
					diff += get_diff(spot_intensities[MASK(p)]) * e;

				if(have_hess)
					hess += e * get_hess(spot_intensities[MASK(p)]) - get_diff(spot_intensities[MASK(p)]).as_col() * get_diff(spot_intensities[MASK(p)]).as_row();
			}
#undef spot_intensities
			log_prob_on = log_prob_on  / (2 * variance) - n_pix * ::log(2*M_PI*variance)/2;
			log_prob_off= log_prob_off / (2 * variance) - n_pix * ::log(2*M_PI*variance)/2;
			diff /= variance;

			log_prob[frame] = make_pair(log_prob_on, log_prob_off);

			if(have_diff)
				diff_log_prob[frame] = diff;

			if(have_hess)
				hess_log_prob[frame] = hess;
		}
	}
		
	double log(int state, int obs) const
	{
		assert(state >=0 && state <= 2);
		assert(obs >=0 && obs < (int)log_prob.size());

		if(state == 0)
			return log_prob[obs].first;
		else
			return log_prob[obs].second;
	}	

	Vector<4>  diff_log(int state, int obs) const
	{
		assert(state >=0 && state <= 2);
		assert(obs >=0 && obs < (int)diff_log_prob.size());

		if(state == 0)
			return diff_log_prob[obs];
		else
			return Zeros;
	}	

	Matrix<4>  hess_log(int state, int obs) const
	{
		assert(state >=0 && state <= 2);
		assert(obs >=0 && obs < (int)hess_log_prob.size());

		if(state == 0)
			return hess_log_prob[obs];
		else
			return Zeros;
	}	

};

#undef SWBG_HAVE_MASK
#undef SWBG_HAVE_DRIFT
#undef SWBG_NAME
#undef SWBG_SPOT_INTENSITIES
#undef MASK

