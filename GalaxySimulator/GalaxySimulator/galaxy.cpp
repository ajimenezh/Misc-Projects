#include "galaxy.h"

#include <random>
#include "math.h"
#include <memory>
#include <functional>
#include <iostream>
#include <chrono>
#include<thread>

#include "Octree.h"
#include "draw_utils.h"

const bool solve_with_threads = true;

double Galaxy::DarkMatterHaloDensity(double r) {
	double q = lambda / r_c;
	static double pi = acos(-1.0);
	double alpha = 1.0 / (1 - sqrt(pi)*q*exp(q*q)*(1 - erf(q)));

	return alpha / (2 * pow(pi, 1.5)*r_c) * exp(-r * r / (r_c*r_c)) / 
		(r*r + lambda * lambda);
}

double Galaxy::StellarBulgeDensity(double r) {
	return r * r / (r + a) / (r + a);
}

void Galaxy::GenerateCoordinates() {

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0, 1.0);

	//for (int i = 0; i < (int) nodes_.size(); i++) {
	//	nodes_[i].x = distribution(generator);
	//	nodes_[i].y = distribution(generator);
	//	nodes_[i].z = distribution(generator);
	//}

	int n = (int)nodes_disk_.size();
	double r = h0;
	for (int i = 0; i < n; i++) {
		double r_n = 10 * distribution(generator) * h0;

		double A = exp(-r_n / h0) / exp(-r / h0);

		double U = distribution(generator);

		if (A >= U) {
			r = r_n;
		}

		double theta = distribution(generator)*2.0*acos(-1.0);
		nodes_disk_[i].x = r * cos(theta);
		nodes_disk_[i].y = r * sin(theta);
	}

	double z = 0.0;
	for (int i = 0; i < n; i++) {
		double z_n = 5 * 2 * (distribution(generator) - 0.5) * z0;

		double A = pow(cosh(z_n / z0), -2) / pow(cosh(z_n / z0), -2);

		double U = distribution(generator);

		if (A >= U) {
			z = z_n;
		}

		nodes_disk_[i].z = z;
	}

	n = (int)nodes_halo_.size();
	r = r_c;
	for (int i = 0; i < n; i++) {
		double r_n = 5 * distribution(generator) * r_c;

		double A = DarkMatterHaloDensity(r_n) / DarkMatterHaloDensity(r);

		double U = distribution(generator);

		if (A >= U) {
			r = r_n;
		}

		double theta = distribution(generator)*2.0*acos(-1.0);
		double phi = distribution(generator)*acos(-1.0);
		nodes_halo_[i].x = r * cos(theta) * sin(phi);
		nodes_halo_[i].y = r * sin(theta) * sin(phi);
		nodes_halo_[i].z = r * cos(phi);
	}

	n = (int)nodes_bulge_.size();
	r = a;
	for (int i = 0; i < n; i++) {
		double r_n = 5 * distribution(generator) * a;

		double A = StellarBulgeDensity(r_n) / StellarBulgeDensity(r);

		double U = distribution(generator);

		if (A >= U) {
			r = r_n;
		}

		double theta = distribution(generator)*2.0*acos(-1.0);
		double phi = distribution(generator)*acos(-1.0);
		nodes_bulge_[i].x = r * cos(theta) * sin(phi);
		nodes_bulge_[i].y = r * sin(theta) * sin(phi);
		nodes_bulge_[i].z = r * cos(phi);
	}
		
}

void Galaxy::GenerateVelocities() {
	GenerateVelocitiesStellarDisk();
	GenerateVelocitiesDarkMatterHalo();
	GenerateVelocitieStellarBulge();
}

void Galaxy::GenerateInitialConditions() {
	GenerateCoordinates();
	//GenerateVelocities();
	ApplyCoordinateTransformations();
}

void Galaxy::ApplyCoordinateTransformations() {

	std::vector<std::vector<double> > rot_x = 
			{ {1.0, 0.0, 0.0}, 
			  {0.0, cos(ang_x), -sin(ang_x)}, 
			  {0.0, sin(ang_x), cos(ang_x)} };

	std::vector<std::vector<double> > rot_y =
			{ {cos(ang_y), 0.0, sin(ang_y)},
			  {0.0, 1.0, 0.0},
			  {-sin(ang_y), 0.0, cos(ang_x)} };

	std::vector<std::vector<double> > rot_z =
			{ {cos(ang_z), -sin(ang_z), 0.0},
			  {sin(ang_z), cos(ang_z), 0.0},
			  {0.0, 0.0, 1.0} };

	auto nodes = { &nodes_disk_ , &nodes_halo_ , &nodes_bulge_ };
	for (auto* it_nodes : nodes) {
		for (auto& node : (*it_nodes)) {
			Vector x(node.x, node.y, node.z);

			x.ApplyRotation(rot_x);
			x.ApplyRotation(rot_y);
			x.ApplyRotation(rot_z);

			node.x = x.x + pos_x;
			node.y = x.y + pos_y;
			node.z = x.z + pos_z;
		}
	}
}

double V(double r, double h0) {
	return 2.0 / (r*h0)*(2 * h0*h0 - exp(-r / h0)*(2 * h0*h0 + 2 * h0*r + r * r));
}

double find_by_interpolation(const std::vector<double>& y, const std::vector<double>& x, double x0) {
	if (x0 <= x[0]) {
		return y[0];
	}
	else if (x0 >= x.back()) {
		return y.back();
	}

	std::vector<double>::const_iterator it = std::lower_bound(x.begin(), x.end(), x0);
	int idx = it - x.begin();

	if (idx == 0) {
		return y[0];
	}

	return (y[idx] - y[idx - 1]) / (x[idx] - x[idx - 1])*(x0 - x[idx - 1]) + y[idx - 1];
}

double calc_derivative(std::function<double(double)> f, double r, double eps) {
	return (f(r + eps) - f(r - eps)) / (2 * eps);
}

double calc_second_derivative(std::function<double(double)> f, double r, double eps) {
	return (f(r + eps) - 2.0*f(r) + f(r - eps)) / (eps * eps);
}

void Galaxy::GenerateVelocitiesStellarDisk() {

	double delta = 1.0e-2;
	double h0 = 1.0;
	double z0 = 1.0;
	double l = 100.0 * h0;

	int n = l / delta;

	double R = l;

	std::vector<double> surf_density(n / 2, 0.0);
	std::vector<double> m_delta(n, 0.0);

	double pi = acos(-1.0);

	for (int i = 0; i < n; i++) {
		double r = i * delta;

		int k = r / delta;
		double sum = 0.0;
		for (int j = 1; j < k; j++) {
			double rr = j * delta;

			sum += V(rr, h0)*rr / sqrt(r*r - rr * rr) * delta;
		}
		m_delta[i] = 2 / pi * sum;
	}

	std::vector<double> radius;
	std::vector<double> vz;
	for (int i = 0; i < n / 2; i++) {
		double r = i * delta;

		int k = (l - r) / delta;
		double sum = 0.0;
		for (int j = 1; j < k; j++) {
			double rr = (i + j) * delta;
			sum += m_delta[j + i] / (rr*rr*sqrt(1 - (r / rr)*(r / rr))) * delta;
		}
		surf_density[i] = 1.0 / 2.0 / pi * sum;
		radius.push_back(r);
		vz.push_back(sqrt(surf_density[i] * z0 * pi));
	}

	//plot(radius, vz);
	//show();

	// V_r(R=h) = 0.6vz(R=0)
	std::vector<double> vr = vz;
	for (int i = 0; i < (int)vr.size(); i++) {
		vr[i] = exp(-radius[i]/h0);
	}
	double factor = vz[0] / find_by_interpolation(vr, radius, h0);
	for (int i = 0; i < (int) vr.size(); i++) {
		vr[i] *= factor;
	}

	std::vector<double> vphi(1);
	std::vector<double> vphi_sigma(1);
	double eps = 1.0e-6;
	for (int i = 1; i < (int)radius.size(); i++) {
		double r = radius[i];
		auto fp = std::bind(&Galaxy::PotentialRAnalyticalR, this, std::placeholders::_1);
		double dphi = calc_derivative(fp, r, 1.0e-6);
		double vc_2 = r * dphi;
		double omega_2 = 1 / r * dphi;
		double vr_2 = find_by_interpolation(vr, radius, r);
		vr_2 = vr_2 * vr_2;
		double k_2 = 3 / r * dphi + calc_second_derivative(fp, r, eps);

		double vr_l = find_by_interpolation(vr, radius, r - eps);
		double vr_r = find_by_interpolation(vr, radius, r + eps);

		double vz_l = find_by_interpolation(vz, radius, r - eps);
		double vz_r = find_by_interpolation(vz, radius, r + eps);

		double arg1 = (log(vr_r*vr_r) - log(vr_l*vr_l)) / (log(r + eps) - log(r - eps));

		double vphi_2 = vc_2 + vr_2 * (1.0 - k_2 / (4 * omega_2) - r/h0 + arg1);

		double sigma = sqrt(vr_2*k_2/(4*omega_2));

		vphi.push_back(sqrt(abs(vphi_2)));
		vphi_sigma.push_back(sigma);
	}

	vphi[0] = vphi[1];
	vphi_sigma[0] = vphi_sigma[1];

	std::default_random_engine generator;
	for (int i = 0; i < (int)nodes_disk_.size(); i++) {

		double r = sqrt(nodes_disk_[i].x*nodes_disk_[i].x +
			nodes_disk_[i].y*nodes_disk_[i].y);

		double vr_dispersion = find_by_interpolation(vr, radius, r);
		std::normal_distribution<double> distribution_r(0.0, vr_dispersion);
		double vr_1 = distribution_r(generator);

		double vz_dispersion = find_by_interpolation(vz, radius, r);
		std::normal_distribution<double> distribution_z(0.0, vz_dispersion);
		double vz_1 = distribution_z(generator);

		double vphi_mean = find_by_interpolation(vphi, radius, r);
		double vphi_dispersion = find_by_interpolation(vphi_sigma, radius, r);
		std::normal_distribution<double> distribution__phi(vphi_mean, vphi_dispersion);
		double vphi_1 = distribution__phi(generator);

		double phi = atan2(nodes_disk_[i].y, nodes_disk_[i].x);
		nodes_disk_[i].vx = vr_1 * cos(phi) - r * vphi_1*sin(phi);
		nodes_disk_[i].vy = vr_1 * sin(phi) - r * vphi_1*cos(phi);
		nodes_disk_[i].vz = vz_1;
	}

}

double DarkMatterHaloVelocityDistribution(double v, double sigma) {
	static double pi = acos(-1.0);
	return 4 * pi*pow(1.0 / (2 * pi*sigma*sigma), 1.5)*v*v*exp(-v * v / (2 * sigma*sigma));
}

void Galaxy::GenerateVelocitiesDarkMatterHalo() {
	double l = r_c * 5;
	int n = 10000;
	double delta = l / n;
	std::vector<double> integrand(n);
	std::vector<double> radius;
	for (int i = 0; i < n; i++) {
		double r = i * delta;
		auto fp = std::bind(&Galaxy::PotentialRAnalyticalR, this, std::placeholders::_1);
		double dphi = calc_derivative(fp, r, 1.0e-6);
		integrand[i] = DarkMatterHaloDensity(r)*dphi;
		radius.push_back(r);
	}

	integrand[n - 1] *= delta;
	for (int i = n - 2; i >= 0; i--) {
		integrand[i] = integrand[i] * delta + integrand[i + 1];
	}

	std::vector<double> vr_2;
	for (int i = 0; i < n; i++) {
		vr_2.push_back(1.0 / DarkMatterHaloDensity(radius[i]) * integrand[i]);
	}

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0, 1.0);

	static double pi = acos(-1.0);

	for (int i = 0; i < (int)nodes_halo_.size(); i++) {
		double r = sqrt(nodes_halo_[i].x*nodes_halo_[i].x +
			nodes_halo_[i].y*nodes_halo_[i].y);

		int n_samples = 100;

		double vr = sqrt(find_by_interpolation(vr_2, radius, r));
		double v = 1.5*sqrt(vr);
		std::vector<double> samples;
		for (int j = 0; j < n_samples; j++) {
			double v_n = 5 * distribution(generator) * vr;

			double A = DarkMatterHaloVelocityDistribution(v_n, vr) / DarkMatterHaloVelocityDistribution(v, vr);

			double U = distribution(generator);

			if (A >= U) {
				v = v_n;
			}

			samples.push_back(v);
		}

		int idx = (int)distribution(generator)*samples.size();
		if (idx == samples.size()) idx--;
		v = samples[idx];
		double phi = distribution(generator) * 2 * pi;
		double theta = distribution(generator) * pi;
		nodes_halo_[i].vx = v * sin(theta)*cos(phi);
		nodes_halo_[i].vy = v * sin(theta)*sin(phi);
		nodes_halo_[i].vz = v * cos(theta);
	}

}

void Galaxy::GenerateVelocitieStellarBulge() {
	double l = r_c * 5;
	int n = 10000;
	double delta = l / n;
	std::vector<double> integrand(n);
	std::vector<double> radius;
	for (int i = 0; i < n; i++) {
		double r = i * delta;
		auto fp = std::bind(&Galaxy::PotentialRAnalyticalR, this, std::placeholders::_1);
		double dphi = calc_derivative(fp, r, 1.0e-6);
		integrand[i] = StellarBulgeDensity(r)*dphi;
		radius.push_back(r);
	}

	integrand[n - 1] *= delta;
	for (int i = n - 2; i >= 0; i--) {
		integrand[i] = integrand[i] * delta + integrand[i + 1];
	}

	std::vector<double> vr_2;
	for (int i = 0; i < n; i++) {
		vr_2.push_back(1.0 / StellarBulgeDensity(radius[i]) * integrand[i]);
	}

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0, 1.0);

	static double pi = acos(-1.0);

	for (int i = 0; i < (int)nodes_bulge_.size(); i++) {
		double r = sqrt(nodes_bulge_[i].x*nodes_bulge_[i].x +
			nodes_bulge_[i].y*nodes_bulge_[i].y);

		int n_samples = 100;

		double vr = sqrt(find_by_interpolation(vr_2, radius, r));
		double v = 1.5*sqrt(vr);
		std::vector<double> samples;
		for (int j = 0; j < n_samples; j++) {
			double v_n = 5 * distribution(generator) * vr;

			double A = DarkMatterHaloVelocityDistribution(v_n, vr) / DarkMatterHaloVelocityDistribution(v, vr);

			double U = distribution(generator);

			if (A >= U) {
				v = v_n;
			}

			samples.push_back(v);
		}

		int idx = (int)distribution(generator)*samples.size();
		if (idx == samples.size()) idx--;
		v = samples[idx];
		double phi = distribution(generator) * 2 * pi;
		double theta = distribution(generator) * pi;
		nodes_bulge_[i].vx = v * sin(theta)*cos(phi);
		nodes_bulge_[i].vy = v * sin(theta)*sin(phi);
		nodes_bulge_[i].vz = v * cos(theta);
	}

}

void Galaxy::PlotGalaxy() {

	std::vector<std::vector<double>> x, y, z;

	int mask = 7;

	if (mask & 1) {
		for (int i = 0; i < (int)nodes_disk_.size(); i++) {
			std::vector<double> x_row, y_row, z_row;

			x_row.push_back(nodes_disk_[i].x);
			y_row.push_back(nodes_disk_[i].y);
			z_row.push_back(nodes_disk_[i].z);

			x.push_back(x_row);
			y.push_back(y_row);
			z.push_back(z_row);
		}
	}

	if (mask & 2) {
		for (int i = 0; i < (int)nodes_halo_.size(); i++) {
			std::vector<double> x_row, y_row, z_row;

			x_row.push_back(nodes_halo_[i].x);
			y_row.push_back(nodes_halo_[i].y);
			z_row.push_back(nodes_halo_[i].z);

			x.push_back(x_row);
			y.push_back(y_row);
			z.push_back(z_row);
		}
	}

	if (mask & 4) {
		for (int i = 0; i < (int)nodes_bulge_.size(); i++) {
			std::vector<double> x_row, y_row, z_row;

			x_row.push_back(nodes_bulge_[i].x);
			y_row.push_back(nodes_bulge_[i].y);
			z_row.push_back(nodes_bulge_[i].z);

			x.push_back(x_row);
			y.push_back(y_row);
			z.push_back(z_row);
		}
	}

	plot(x, y, z, -10, 10);
	show();
}

double SquareDistance(const Node& a, const Node& b) {
	return (a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y) + (a.z - b.z)*(a.z - b.z);
}

/*
* The potential for a particle i is the cotribution of each particle.
*/
double Galaxy::Potential(const Node& nd) {
	double phi = 0;
	double eps = 1.0;
	double m = 1.0 / (int)nodes_disk_.size();
	for (int i = 0; i < (int)nodes_disk_.size(); i++) {
		phi -= m / sqrt(SquareDistance(nodes_disk_[i], nd) + eps * eps);
	}

	m = 1.0 / (int)nodes_halo_.size();
	for (int i = 0; i < (int)nodes_halo_.size(); i++) {
		phi -= m / sqrt(SquareDistance(nodes_halo_[i], nd) + eps * eps);
	}

	m = 1.0 / (int)nodes_bulge_.size();
	for (int i = 0; i < (int)nodes_bulge_.size(); i++) {
		phi -= m / sqrt(SquareDistance(nodes_bulge_[i], nd) + eps * eps);
	}

	return phi;
}

double Galaxy::PotentialRAnalytical(const Node& nd) {
	double phi = 0;
	double eps = 1.0;
	double pi = acos(-1.0);
	
	double r = sqrt(nd.x*nd.x + nd.y*nd.y);

	double q = lambda / r_c;
	double alpha = 1.0 / (1 - sqrt(pi)*q*exp(q*q)*(1 - erf(q)));

	// TODO: The factor 1.0/r is wrong, it should be M(r)/r.
	phi = 2 * exp(-r / h0) - 1.0/(r + a) - 1.0/r + 
		alpha/(sqrt(pi)*r_c)*exp(-(r/r_c)*(r / r_c) - q*q);

	return phi;
}

void Galaxy::PlotPotential() {

	std::vector<double> r, phi, phi2;

	int n = 1000;
	double max_r = 10.0 * h0;
	double delta = max_r / n;

	for (int i = 1; i <= n; i++) {
		double rr = i * delta;
		Node nd;
		nd.x = rr;

		r.push_back(rr);
		phi.push_back(Potential(nd));
	}

	for (int i = 1; i <= n; i++) {
		double rr = i * delta;
		Node nd;
		nd.x = rr;

		phi2.push_back(PotentialRAnalytical(nd));
	}
		 
	plot(r, phi);
	plot(r, phi2);
	show();
}

std::vector<tIntegrationMethods> integrationScheme1 = { 
	HALF_STEP_POSITION_INTEGRATION, 
	STEP_VELOCITY_INTEGRATION,
	HALF_STEP_POSITION_INTEGRATION 
};

std::vector<tIntegrationMethods> integrationScheme2 = {
	HALF_STEP_VELOCITY_INTEGRATION,
	STEP_POSITION_INTEGRATION,
	HALF_STEP_VELOCITY_INTEGRATION
};

void Universe::PartialSolve(int s, int e, int n, double delta_t, Octree* octree, tIntegrationMethods integration_method) {
	for (int i = s; i < e; i++) {
		Vector f = -octree->CalcForce(nodes_[i]);

		switch (integration_method) {
			case HALF_STEP_POSITION_INTEGRATION:
				(*nodes_[i]).x += (*nodes_[i]).vx*delta_t / 2;
				(*nodes_[i]).y += (*nodes_[i]).vy*delta_t / 2;
				(*nodes_[i]).z += (*nodes_[i]).vz*delta_t / 2;
				break;
			case STEP_VELOCITY_INTEGRATION:
				(*nodes_[i]).vx += f.x * delta_t;
				(*nodes_[i]).vy += f.y * delta_t;
				(*nodes_[i]).vz += f.z * delta_t;
				break;
			case STEP_POSITION_INTEGRATION:
				(*nodes_[i]).x += (*nodes_[i]).vx*delta_t;
				(*nodes_[i]).y += (*nodes_[i]).vy*delta_t;
				(*nodes_[i]).z += (*nodes_[i]).vz*delta_t;
				break;
			case HALF_STEP_VELOCITY_INTEGRATION:
				(*nodes_[i]).vx += f.x * delta_t / 2;
				(*nodes_[i]).vy += f.y * delta_t / 2;
				(*nodes_[i]).vz += f.z * delta_t / 2;
				break;
		}
	}
}

void Universe::Solve(int n, double delta_t) {
	auto start = std::chrono::high_resolution_clock::now();
	int mod = 10;
	int numthreads = std::thread::hardware_concurrency();

	for (int it = 0; it < n; it++) {
		if (it!= 0 && it % mod == 0) {
			auto finish = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = finish - start;

			std::cout << "Iteration: " << it << std::endl;
			std::cout << "Elapsed time: " << elapsed.count() << " s\n";
			start = std::chrono::high_resolution_clock::now();
		}

		Octree octree;
		octree_cache.FreeOctree();
		octree.BuildTree({ NodesBegin(), NodesEnd() });

		if (solve_with_threads) {
			size_t start = 0;
			size_t end = nodes_.size();
			size_t chunk = (numthreads > 0 ? (end - start + (numthreads - 1)) / numthreads : end - start);

			std::vector<std::thread> t(numthreads);
			for (int k = 0; k < integrationScheme1.size(); k++) {
				for (int i = 0; i < numthreads; i++)
				{
					int s = start + i * chunk;
					int e = std::min(s + chunk, end);
					t[i] = std::thread(&Universe::PartialSolve, this, s, e, n, delta_t, &octree, integrationScheme1[k]);
				}
				for (int i = 0; i < numthreads; i++)
					t[i].join();
			}
		}
		else {
			for (auto* node : nodes_) {
				//Vector f2 = Force(*node);

				Vector f = -octree.CalcForce(node);

				//f = f / (*nodes)[i].m*delta_t;
				f = f * delta_t;
				(*node).vx += f.x;
				(*node).vy += f.y;
				(*node).vz += f.z;

				(*node).x += (*node).vx*delta_t;
				(*node).y += (*node).vy*delta_t;
				(*node).z += (*node).vz*delta_t;
			}
		}
	}
}

Vector Universe::Force(const Node& nd) {

	Vector f;
	for (auto* node : nodes_) {
		if (*node != nd) {
			Vector v(nd.x - (*node).x, nd.y - (*node).y, nd.z - (*node).z);
			double tmp = v.sqr() + eps * eps;
			f += v * (*node).m / (tmp*sqrt(tmp));
		}
	}

	return f;
}

void Universe::Plot() {

	std::vector<std::vector<double>> x, y, z;

	for (int i = 0; i < (int)nodes_.size(); i++) {
		std::vector<double> x_row, y_row, z_row;

		x_row.push_back((*nodes_[i]).x);
		y_row.push_back((*nodes_[i]).y);
		z_row.push_back((*nodes_[i]).z);

		x.push_back(x_row);
		y.push_back(y_row);
		z.push_back(z_row);
	}

	plot(x, y, z, -10, 50);
	show();
}