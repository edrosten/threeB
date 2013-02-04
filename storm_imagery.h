#ifndef STORM_INCLUDE_STORM_IMAGERY_H
#define STORM_INCLUDE_STORM_IMAGERY_H
#include <cvd/image.h>
#include <cvd/byte.h>
#include <utility>
#include <vector>

//! @cond Doxygen_Suppress
std::vector<CVD::Image<float> > load_and_preprocess_images2(const std::vector<std::string>& names);
std::vector<CVD::Image<float> > load_and_preprocess_images(const std::vector<std::string>& names);
CVD::Image<float> preprocess_image(const CVD::Image<float>&);
std::pair<float, float> mean_and_variance(const std::vector<CVD::Image<float> >& images);
std::vector<CVD::Image<float> > load_and_normalize_images(const std::vector<std::string>& files);
//! @endcond
#endif



