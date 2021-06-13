#pragma once

const double eps = 1.0e-8;

struct Node {
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;

	double vx = 0.0;
	double vy = 0.0;
	double vz = 0.0;

	double m = 0.0;

	bool operator!=(const Node& nd) const {
		if (abs(nd.x - x) > eps) return true;
		if (abs(nd.y - y) > eps) return true;
		if (abs(nd.z - z) > eps) return true;
		return false;
	}
};

struct Vector {
	Vector() {}
	Vector(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

	double sqr() {
		return x * x + y * y + z * z;
	}

	void ApplyRotation(const std::vector<std::vector<double> >& rot) {
		std::vector<double> v = { x, y, z };
		std::vector<double> tmp(3, 0.0);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				tmp[i] += rot[i][j] * v[j];
			}
		}
		x = tmp[0];
		y = tmp[1];
		z = tmp[2];
	}

	void operator+=(const Vector& v) {
		x += v.x;
		y += v.y;
		z += v.z;
	}

	void operator-=(const Vector& v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;
	}

	Vector& operator*(double val) {
		x *= val;
		y *= val;
		z *= val;
		return *this;
	}

	Vector& operator/(double val) {
		x /= val;
		y /= val;
		z /= val;
		return *this;
	}

	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
};

inline Vector operator-(const Vector& v) {
	Vector tmp = v;
	tmp.x = -tmp.x;
	tmp.y = -tmp.y;
	tmp.z = -tmp.z;
	return tmp;
}