#pragma once
#include <vector>
#include "utils.h"

const int MAX_NODES = 200000;
const int NODES_PER_LEAF = 8;
const double inf = 1 << 20;

struct ListData  {
	std::vector<Node*>::iterator list_begin;
	std::vector<Node*>::iterator list_end;

	std::pair<double, double> x = {inf, -inf};
	std::pair<double, double> y = { inf, -inf };
	std::pair<double, double> z = { inf, -inf };

	double m = 0.0;

	double mx = 0.0;
	double my = 0.0;
	double mz = 0.0;

	void clear() {
		x = { inf, -inf };
		y = { inf, -inf };
		z = { inf, -inf };

		m = 0.0;

		mx = 0.0;
		my = 0.0;
		mz = 0.0;
	}
};

class Octree {
public:
	Octree();
	~Octree();

	void BuildTree(ListData data);

	std::vector<ListData> DivideByX(std::vector<ListData> lists);
	std::vector<ListData> DivideByY(std::vector<ListData> lists);
	std::vector<ListData> DivideByZ(std::vector<ListData> lists);

	void CalcData();

	Vector CalcForce(Node* nd);

	void Clear() {
		for (int i = 0; i < 8; i++) children[i] = NULL;
		parent = NULL;
		objects.clear();
	}

private:
	Octree* children[8];
	Octree* parent;
	ListData objects;
};

bool CheckCondition(Node* nd, ListData& data);

class OctreeCache {
public:
	OctreeCache() {
		n = MAX_NODES;
		octree_nodes.resize(n);
		for (int i = 0; i < n; i++) {
			octree_nodes[i] = new Octree();
		}
		free_idx = 0;
	}
	~OctreeCache() {
		for (int i = 0; i < n; i++) {
		delete octree_nodes[i];
		}
	}
	Octree* GetOctree() {
		if (free_idx == n) throw("Octree nodes limit");
		Octree* x = octree_nodes[free_idx++];
		x->Clear();
		return x;
	}
	void FreeOctree(Octree* node) {
		std::swap(octree_nodes[free_idx - 1], node);
		free_idx--;
	}
	void FreeOctree() {
		free_idx = 0;
	}
private:
	std::vector<Octree*> octree_nodes;
	int n;
	int free_idx;
};

static OctreeCache octree_cache;
