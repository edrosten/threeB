#ifndef MULTISPOT5_PLACE_METHODS_H
#define MULTISPOT5_PLACE_METHODS_H
#include "multispot5.h"

StateParameters place_spots_uniform(int num_spots, const std::vector<CVD::ImageRef>& pixels, const CVD::ImageRef& size);
StateParameters place_spots_intensity_sampled(int num_spots, const std::vector<CVD::ImageRef>& pixels, const std::vector<CVD::Image<float> >&);
#endif
