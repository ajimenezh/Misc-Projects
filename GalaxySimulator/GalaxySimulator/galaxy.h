#pragma once
#include "string.h"
#include <vector>
#include "utils.h"

class Octree;

enum tIntegrationMethods {
	HALF_STEP_POSITION_INTEGRATION,
	STEP_POSITION_INTEGRATION,
	HALF_STEP_VELOCITY_INTEGRATION,
	STEP_VELOCITY_INTEGRATION
};

class Galaxy {
public:

	Galaxy() {
	}

	void SetName(const std::string& name) {
		name_ = name;
	}

	void SetCenter(double x, double y, double z) {
		pos_x = x;
		pos_y = y;
		pos_z = z;
	}

	void SetAngles(double x, double y, double z) {
		ang_x = x;
		ang_y = y;
		ang_z = z;
	}

	void SetMass(double x, double y, double z) {
		mass_disk = x;
		mass_halo = y;
		mass_bulge = z;
	}

	void SetMassParticles() {
		for (auto& node : nodes_disk_) {
			node.m = mass_disk / nodes_disk_.size();
		}
		for (auto& node : nodes_halo_) {
			node.m = mass_halo / nodes_halo_.size();
		}
		for (auto& node : nodes_bulge_) {
			node.m = mass_bulge / nodes_bulge_.size();
		}
	}

	void SetNumParticlesDisk(int n) {
		nodes_disk_.resize(n);
	}

	void SetNumParticlesHalo(int n) {
		nodes_halo_.resize(n);
	}

	void SetNumParticlesBulge(int n) {
		nodes_bulge_.resize(n);
	}

	void SetH0(double x) {
		h0 = x;
	}

	void SetZ0(double x) {
		z0 = x;
	}

	void SetLambda(double x) {
		lambda = x;
	}

	void SetRc(double x) {
		r_c = x;
	}

	void SetA(double x) {
		a = x;
	}

	void GenerateInitialConditions();

	void GenerateCoordinates();

	void GenerateVelocities();

	void ApplyCoordinateTransformations();

	void GenerateVelocitiesStellarDisk();

	void GenerateVelocitiesDarkMatterHalo();

	void GenerateVelocitieStellarBulge();

	void PlotGalaxy();

	void PlotPotential();

	double Potential(const Node& nd);

	double PotentialRAnalyticalR(double r) {
		Node nd;
		nd.x = r;
		return PotentialRAnalytical(nd);
	}

	double PotentialRAnalytical(const Node& nd);

	void BuildNodes(std::vector<Node*>& nodes) {
		for (auto& node : nodes_disk_) {
			nodes.push_back(&node);
		}
		for (auto& node : nodes_halo_) {
			nodes.push_back(&node);
		}
		for (auto& node : nodes_bulge_) {
			nodes.push_back(&node);
		}
	}

private:

	double DarkMatterHaloDensity(double r);

	double StellarBulgeDensity(double r);

private:
	std::string name_;
	std::vector<Node> nodes_disk_;
	std::vector<Node> nodes_halo_;
	std::vector<Node> nodes_bulge_;

	double mass_disk;
	double mass_halo;
	double mass_bulge;

	double h0;
	double z0;
	double lambda;
	double r_c;
	double a;

	double pos_x;
	double pos_y;
	double pos_z;

	double ang_x;
	double ang_y;
	double ang_z;
};

class Universe {
public:

	Universe(int n) {
		galaxies.resize(n + 1);
		size = n;
	}

	Galaxy& operator[](int i) {
		return galaxies[i];
	}

	typedef Galaxy* iterator;
	typedef const Galaxy* const_iterator;

	iterator begin() { return &galaxies[0]; }
	const_iterator begin() const { return &galaxies[0]; }
	iterator end() { return &galaxies[size]; }
	const_iterator end() const { return &galaxies[size]; }

	void BuildNodes() {
		for (auto& galaxy : galaxies) {
			galaxy.BuildNodes(nodes_);
		}
	}

	std::vector<Node*>::iterator NodesBegin() {
		if (nodes_.size() == 0) {
			BuildNodes();
		}
		return nodes_.begin();
	}

	std::vector<Node*>::iterator NodesEnd() {
		if (nodes_.size() == 0) {
			BuildNodes();
		}
		return nodes_.end();
	}

	void Plot();

	void Solve(int n, double delta_t);

	void PartialSolve(int s, int e, int n, double delta_t, Octree* octree, tIntegrationMethods integration_method);

	Vector Force(const Node& nd);

private:
	std::vector<Galaxy> galaxies;
	std::size_t size;

	std::vector<Node*> nodes_;
};