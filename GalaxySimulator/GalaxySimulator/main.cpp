#include <iostream>
#include "config_file_reader.h"
#include "draw_utils.h"
#include "Octree.h"
#include "galaxy_builder.h"

void SolveGalaxy() {
	auto config_data = ReadConfigFile("config.txt");

	Universe universe = GalaxyBuilder(config_data);

	for (auto& galaxy : universe) {
		galaxy.GenerateInitialConditions();
		// galaxy.Solve(100, 1.0e-3);
	}

	universe.Solve(100, 5.0e-3);

	universe.Plot();
	// galaxies[0].PlotPotential();
}

double V2(double r, double h0) {
	return 2.0 / (r*h0)*(2*h0*h0 - exp(-r / h0)*(2 * h0*h0 + 2 * h0*r + r * r));
}

void test() {
	double delta = 1.0e-2;
	double h0 = 1.0;
	double z0 = 1.0;
	double l = 100.0 * h0;

	int n = l / delta;

	double R = l;

	std::vector<double> surf_density(n/2, 0.0);
	std::vector<double> m_delta(n, 0.0);

	double pi = acos(-1.0);

	for (int i = 0; i < n; i++) {
		double r = i * delta;

		int k = r / delta;
		double sum = 0.0;
		for (int j = 1; j < k; j++) {
			double rr = j * delta;

			sum += V2(rr, h0)*rr / sqrt(r*r - rr*rr);
		} 
		m_delta[i] = 2 / pi * sum;
	}

	std::vector<double> y;
	std::vector<double> vz;
	for (int i = 0; i < n/2; i++) {
		double r = i * delta;

		int k = (l - r) / delta;
		double sum = 0.0;
		for (int j = 1; j < k; j++) {
			double rr = (i+j) * delta;
			sum += m_delta[j + i] / (rr*rr*sqrt(1 - (r / rr)*(r / rr)));
		}
		surf_density[i] = 1.0 / 2.0 / pi * sum;
		y.push_back(r);
		vz.push_back(sqrt(surf_density[i]*z0*pi));
	}

	plot(y, vz);
	show();
}

int main() {

	SolveGalaxy();

	// test();

	return 0;
}