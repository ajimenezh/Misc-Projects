#include "draw_utils.h"
#include "matplotlibcpp.h"

void plot(const std::vector<double>& x, const std::vector<double> y) {
	matplotlibcpp::plot(x, y);
}

void plot(const std::vector<std::vector<double>>& x, const std::vector<std::vector<double>>& y, 
	const std::vector<std::vector<double>>& z, double lim_min, double lim_max) {
	std::map<std::string, std::string> params;

	params["marker"] = ".";

	matplotlibcpp::scatter_3d(x, y, z, params);
	matplotlibcpp::xlim(lim_min, lim_max);
	matplotlibcpp::ylim(lim_min, lim_max);
}

void show() {
	matplotlibcpp::show();
}