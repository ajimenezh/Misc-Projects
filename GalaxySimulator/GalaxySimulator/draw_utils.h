#pragma once
#include <vector>

void plot(const std::vector<double>& x, const std::vector<double> y);

void plot(const std::vector<std::vector<double>>& x, const std::vector<std::vector<double>>& y, 
	const std::vector<std::vector<double>>& z, double lim_min, double lim_max);

void show();