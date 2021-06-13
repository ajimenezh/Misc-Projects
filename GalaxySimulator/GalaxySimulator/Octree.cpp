#include "Octree.h"
#include <algorithm>


Octree::Octree() {
	for (int i = 0; i < 8; i++) children[i] = NULL;
	parent = NULL;
} 


Octree::~Octree()
{
}

void Octree::BuildTree(ListData data) {
	objects = data;
	if (objects.list_end - objects.list_begin <= NODES_PER_LEAF) {
		CalcData();
		return;
	}

	std::vector<ListData> tmp_objects = DivideByX({ objects });
	tmp_objects = DivideByY({ tmp_objects });
	tmp_objects = DivideByZ({ tmp_objects });

	for (int i = 0; i < 8; i++) {
		Octree* node = octree_cache.GetOctree();
		node->BuildTree(tmp_objects[i]);
		children[i] = node;
	}
	CalcData();
}

std::vector<ListData> Octree::DivideByX(std::vector<ListData> lists) {
	std::vector<ListData> new_objects;
	for (const auto& list : lists) {
		std::sort(list.list_begin, list.list_end,
			[](const Node* a, const Node* b)
		{
			return a->x < b->x;
		});

		int n = (int)(list.list_end - list.list_begin);

		new_objects.emplace_back(ListData({ list.list_begin, list.list_begin + n / 2 }));
		new_objects.emplace_back(ListData({ list.list_begin + n / 2, list.list_end }));
	}

	return new_objects;
}

std::vector<ListData> Octree::DivideByY(std::vector<ListData> lists) {
	std::vector<ListData> new_objects;
	for (const auto& list : lists) {
		std::sort(list.list_begin, list.list_end,
			[](const Node* a, const Node* b)
		{
			return a->y < b->y;
		});

		int n = (int)(list.list_end - list.list_begin);

		new_objects.emplace_back(ListData({ list.list_begin, list.list_begin + n / 2 }));
		new_objects.emplace_back(ListData({ list.list_begin + n / 2, list.list_end }));
	}

	return new_objects;
}

std::vector<ListData> Octree::DivideByZ(std::vector<ListData> lists) {
	std::vector<ListData> new_objects;
	for (const auto& list : lists) {
		std::sort(list.list_begin, list.list_end,
			[](const Node* a, const Node* b)
		{
			return a->z < b->z;
		});

		int n = (int)(list.list_end - list.list_begin);

		new_objects.emplace_back(ListData({ list.list_begin, list.list_begin + n / 2 }));
		new_objects.emplace_back(ListData({ list.list_begin + n / 2, list.list_end }));
	}

	return new_objects;
}

void Octree::CalcData() {
	if (objects.list_begin == objects.list_end) {
		return;
	}

	if (objects.list_end - objects.list_begin <= NODES_PER_LEAF) {
		for (auto it = objects.list_begin; it != objects.list_end; it++) {
			objects.x.first = std::min(objects.x.first, (*it)->x);
			objects.x.second = std::max(objects.x.second, (*it)->x);

			objects.y.first = std::min(objects.y.first, (*it)->y);
			objects.y.second = std::max(objects.y.second, (*it)->y);

			objects.z.first = std::min(objects.z.first, (*it)->z);
			objects.z.second = std::max(objects.z.second, (*it)->z);

			objects.m += (*it)->m;

			objects.mx += (*it)->x*(*it)->m;
			objects.my += (*it)->y*(*it)->m;
			objects.mz += (*it)->z*(*it)->m;
		}

		objects.mx /= objects.m;
		objects.my /= objects.m;
		objects.mz /= objects.m;

		return;
	}

	for (int i = 0; i < 8; i++) {
		if (children[i] != NULL) {
			objects.x.first = std::min(objects.x.first, children[i]->objects.x.first);
			objects.x.second = std::max(objects.x.second, children[i]->objects.x.second);

			objects.y.first = std::min(objects.y.first, children[i]->objects.y.first);
			objects.y.second = std::max(objects.y.second, children[i]->objects.y.second);

			objects.z.first = std::min(objects.z.first, children[i]->objects.z.first);
			objects.z.second = std::max(objects.z.second, children[i]->objects.z.second);

			objects.m += children[i]->objects.m;

			objects.mx += children[i]->objects.mx*children[i]->objects.m;
			objects.my += children[i]->objects.my*children[i]->objects.m;
			objects.mz += children[i]->objects.mz*children[i]->objects.m;
		}
	}

	objects.mx /= objects.m;
	objects.my /= objects.m;
	objects.mz /= objects.m;
}

Vector Octree::CalcForce(Node* nd) {
	if ((objects.list_end - objects.list_begin) <= NODES_PER_LEAF) {
		Vector F;
		for (auto it = objects.list_begin; it != objects.list_end; it++) {
			if (*(*it) != *nd) {
				Vector v(nd->x - (*it)->x, nd->y - (*it)->y, nd->z - (*it)->z);
				double tmp = v.sqr() + eps * eps;
				F += v * (*it)->m / (tmp*sqrt(tmp));
			}
		}
		return F;
	}

	Vector F;
	for (int i = 0; i < 8; i++) {
		if (children[i] != NULL) {
			if (CheckCondition(nd, children[i]->objects)) {
				Vector v(nd->x - children[i]->objects.mx, 
					nd->y - children[i]->objects.my,
					nd->z - children[i]->objects.mz);
				double tmp = v.sqr() + eps * eps;
				F += v * children[i]->objects.m / (tmp*sqrt(tmp));
			}
			else {
				F += children[i]->CalcForce(nd);
			}
		}
	}

	return F;
}

bool CheckCondition(Node* nd, ListData& data) {
	double d = std::max(std::max(data.x.second - data.x.first,
		data.y.second - data.y.first), data.z.second - data.z.first);
	double r = sqrt((nd->x - data.mx)*(nd->x - data.mx) +
		(nd->y - data.my)*(nd->y - data.my) + (nd->z - data.mz)*(nd->z - data.mz));

	static double theta = 1.0;

	return d / r < theta;
}