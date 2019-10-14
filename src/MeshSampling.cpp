#include "stdafx.h"
#include "MeshSampling.h"
#include <array>

MeshSampling::MeshSampling()
{
}


MeshSampling::~MeshSampling()
{
}

//assemble cc bbx
void MeshSampling::assemble_bbx() {
	auto update_bbx = [&](const double *pos, BlackMesh::BlackMesh<double>::BBX &bbx) {
		if (bbx.xmin > pos[0]) {
			bbx.xmin = pos[0];
		}
		if (bbx.ymin > pos[1]) {
			bbx.ymin = pos[1];
		}
		if (bbx.zmin > pos[2]) {
			bbx.zmin = pos[2];
		}
		if (bbx.xmax < pos[0]) {
			bbx.xmax = pos[0];
		}
		if (bbx.ymax < pos[1]) {
			bbx.ymax = pos[1];
		}
		if (bbx.zmax < pos[2]) {
			bbx.zmax = pos[2];
		}
	};

	cc_bbx.clear();
	cc_bbx.resize(msh->GetNumComponents(), BlackMesh::BlackMesh<double>::BBX());

	for (int i = 0; i < cc_triid.size(); i++) {
		for (int j = 0; j < cc_triid[i].size(); j++) {
			auto tri = msh->GetTriangles()[cc_triid[i][j]];
			update_bbx(&(msh->GetVertices()[tri.vertices[0]].pos[0]), (cc_bbx[i]));
			update_bbx(&(msh->GetVertices()[tri.vertices[1]].pos[0]), (cc_bbx[i]));
			update_bbx(&(msh->GetVertices()[tri.vertices[2]].pos[0]), (cc_bbx[i]));
		}
	}
}



//assemble data
void MeshSampling::assemble_data() {
	auto getTriangleArea = [](const BlackMesh::BlackMesh<double> *msh, const int &id) {
#define x 0
#define y 1
#define z 2
		auto tri = msh->GetTriangles()[id];
		auto p0 = msh->GetVertices()[tri.vertices[0]].pos;
		auto p1 = msh->GetVertices()[tri.vertices[1]].pos;
		auto p2 = msh->GetVertices()[tri.vertices[2]].pos;

		return 0.5*abs(((p1[y] - p0[y])*(p2[z] - p0[z]) + (p1[z] - p0[z])*(p2[x] - p0[x]) + (p1[x] - p0[x])*(p2[y] - p0[y])) -
			((p2[y] - p0[y])*(p1[z] - p0[z]) + (p2[z] - p0[z])*(p1[x] - p0[x]) + (p2[x] - p0[x])*(p1[y] - p0[y])));
#undef x
#undef y
#undef z
	};

	cc_total_area.clear();
	cc_total_area.resize(msh->GetNumComponents(), 0);

	cc_triid.clear();
	cc_triid.resize(msh->GetNumComponents(), vector<int>{});

	face_area.clear();
	face_area.reserve(msh->GetNumTriangles());

	total_area = 0;

	for (int i = 0; i < msh->GetNumTriangles(); i++) {
		double ara = getTriangleArea(msh, i);
		auto tri = msh->GetTriangles()[i];

		cc_total_area[tri.component_id] += ara;
		face_area.push_back(ara);

		cc_triid[tri.component_id].push_back(i);

		total_area += ara;
	}

}

//generate rand points inside a triangle
void MeshSampling::generate_n_rand_seed(BlackMesh::BlackMesh<double> *msh, int id, int n, std::vector<std::array<double, 3>> &outputlist) {
	auto tri = msh->GetTriangles()[id];
	auto p0 = msh->GetVertices()[tri.vertices[0]].pos;
	auto p1 = msh->GetVertices()[tri.vertices[1]].pos;
	auto p2 = msh->GetVertices()[tri.vertices[2]].pos;

	outputlist.clear();
	outputlist.reserve(n);
	for (int i = 0; i < n; i++) {
		double x = (double)rand() / (RAND_MAX + 1.0);
		double y = (double)rand() / (RAND_MAX + 1.0);

		double cof1 = 1 - sqrt(x);
		double cof2 = sqrt(x)*(1 - y);
		double cof3 = y*sqrt(x);

		outputlist.push_back(array<double, 3>{cof1*p0[0] + cof2*p1[0] + cof3*p2[0],
			cof1*p0[1] + cof2*p1[1] + cof3*p2[1],
			cof1*p0[2] + cof2*p1[2] + cof3*p2[2]});

	}
}

//sample  P=max{A/(diag*diag)*N, Pfix}
void MeshSampling::generate_sampling_points(const int total_points, const int minimal_points) {//total points is for all cc
	if (msh->GetNumComponents() == 0) {
		msh->mark_component_with_coherence();
	}
	assemble_data();
	assemble_bbx();




	double diag_len = msh->bbx.update_diag_dist();
	cc_sampled_triid.clear();
	cc_sampled_triid.resize(msh->GetNumComponents(), std::vector<std::pair<int, array<double, 3>>>{});

	for (int i = 0; i < cc_triid.size(); i++) {

		int cc_total = max((int)(cc_total_area[i] / (diag_len*diag_len)*total_points), minimal_points);

		double carry_area = 0;

		for (int j = 0; j < cc_triid[i].size(); j++) {
			
			int triidx = cc_triid[i][j];
			double raw_data = face_area[triidx] / cc_total_area[i] * cc_total + carry_area + 1e-6;
			int pts = (int)raw_data;
			carry_area = raw_data - (double)pts;


			if (pts > 0) {
				std::vector<array<double, 3>> list;
				generate_n_rand_seed(msh, triidx, pts, list);
				for (int mm = 0; mm < list.size(); mm++)
					cc_sampled_triid[i].push_back(std::make_pair(triidx, list[mm]));
			}
		}


	}

}
