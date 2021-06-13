#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <variant>

const int MAX_NUMBER_GALAXIES = 100;

class CustomVariant {
public:
	CustomVariant(double d_) {
		d = d_;
		type = DOUBLE;
	}

	CustomVariant(std::string s_) {
		s = s_;
		type = STRING;
	}

	std::string s;
	double d;

	enum t_type {
		DOUBLE,
		STRING
	};

	t_type type;
};

template<class T>
T* get_if(CustomVariant* obj) {
	std::string tname = typeid(T).name();

	if (tname == "double") {
		return (T*)&obj->d;
	}
	else {
		return (T*)&obj->s;
	}
}

typedef std::vector<std::map<std::string, std::vector<CustomVariant > > > tConfigData;

tConfigData ReadConfigFile(const char* file_name) {
	std::ifstream infile;

	tConfigData data;

	infile.open(file_name);
	if (!infile) {
		std::cout << "Error reading config file" << std::endl;
		return {};
	}

	std::string line;
	int idx = -1;
	while (std::getline(infile, line)) {
		if (line.size() == 0 || line[0] == '#') {
			continue;
		}
		std::istringstream iss(line);
		
		std::string name;
		iss >> name;

		if (name == "new_galaxy") {
			idx++;
			data.resize(idx + 1);
		}
		else if (idx != -1) {
			std::string s;
			while (iss >> s) {
				if (s.length() > 0 && isdigit(s[0])) {
					data[idx][name].push_back(stof(s));
				}
				else {
					data[idx][name].push_back(s);
				}
			}
		}
	}

	infile.close();

	return data;
}
