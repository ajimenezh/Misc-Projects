#pragma once
#include <vector>
#include <map>
#include <variant>
#include "galaxy.h"

Universe GalaxyBuilder(
	std::vector<std::map<std::string, std::vector<CustomVariant > > >& config_data) {

	Universe galaxies(config_data.size());

	for (int i = 0; i < config_data.size(); i++) {
		{
			auto pval = get_if<std::string>(&config_data[i]["name"][0]);
			if (pval) {
				galaxies[i].SetName(*pval);
			}
		}

		{
			auto pval = get_if<double>(&config_data[i]["num_particles_disk"][0]);
			if (pval) {
				galaxies[i].SetNumParticlesDisk(*pval);
			}
		}

		{
			auto pval = get_if<double>(&config_data[i]["num_particles_halo"][0]);
			if (pval) {
				galaxies[i].SetNumParticlesHalo(*pval);
			}
		}

		{
			auto pval = get_if<double>(&config_data[i]["num_particles_bulge"][0]);
			if (pval) {
				galaxies[i].SetNumParticlesBulge(*pval);
			}
		}

		{
			auto pval = get_if<double>(&config_data[i]["h0"][0]);
			if (pval) {
				galaxies[i].SetH0(*pval);
			}
		}

		{
			auto pval = get_if<double>(&config_data[i]["z0"][0]);
			if (pval) {
				galaxies[i].SetZ0(*pval);
			}
		}

		{
			auto pval = get_if<double>(&config_data[i]["lambda"][0]);
			if (pval) {
				galaxies[i].SetLambda(*pval);
			}
		}

		{
			auto pval = get_if<double>(&config_data[i]["r_c"][0]);
			if (pval) {
				galaxies[i].SetRc(*pval);
			}
		}

		{
			auto pval = get_if<double>(&config_data[i]["a"][0]);
			if (pval) {
				galaxies[i].SetA(*pval);
			}
		}

		{
			auto pval = get_if<double>(&config_data[i]["mass"][0]);
			std::vector<double> pos;
			for (int j = 0; j < 3; j++) {
				auto pval = get_if<double>(&config_data[i]["mass"][j]);
				pos.push_back(*pval);
			}
			galaxies[i].SetMass(pos[0], pos[1], pos[2]);
			galaxies[i].SetMassParticles();
		}

		{
			auto pval = get_if<double>(&config_data[i]["position"][0]);
			std::vector<double> pos;
			for (int j = 0; j < 3; j++) {
				auto pval = get_if<double>(&config_data[i]["position"][j]);
				pos.push_back(*pval);
			}
			galaxies[i].SetCenter(pos[0], pos[1], pos[2]);
		}

		{
			auto pval = get_if<double>(&config_data[i]["angles"][0]);
			std::vector<double> pos;
			for (int j = 0; j < 3; j++) {
				auto pval = get_if<double>(&config_data[i]["angles"][j]);
				pos.push_back(*pval);
			}
			galaxies[i].SetAngles(pos[0], pos[1], pos[2]);
		}
	}

	return galaxies;
}