#ifndef MULTISPOT5_PLACE_CHOICE_H
#define MULTISPOT5_PLACE_CHOICE_H
#include "multispot5.h"

void place_and_fit_spots(const std::vector<CVD::Image<float> >& ims, const std::vector<CVD::ImageRef>& region, const CVD::Image<double>& log_ratios, std::string save_spots_file, FitSpotsGraphics& g, const std::string&s="");
void place_and_fit_spots(const std::vector<CVD::Image<float> >& ims, const std::vector<CVD::ImageRef>& region, const CVD::Image<double>& log_ratios, std::ofstream& save_spots_file, FitSpotsGraphics& g, UserInterfaceCallback&, const std::string&s="");

#endif
