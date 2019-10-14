#pragma once
#include <array>
#include "BlackMesh.h"

class MeshSampling
{
private:
	BlackMesh::BlackMesh<double> *msh;


	//generate rand points inside a triangle
	void generate_n_rand_seed(BlackMesh::BlackMesh<double> *msh, int id, int n,
		std::vector<std::array<double, 3>> &outputlist);

	//assemble cc bbx
	void assemble_bbx();

	//assemble data
	void assemble_data();

public:

	std::vector<std::vector<std::pair<int, array<double, 3>>>> cc_sampled_triid;//ccid-triangle id, pos

	std::vector<double> cc_total_area;
	std::vector<vector<int>> cc_triid;
	std::vector<double> face_area;

	double total_area;

	std::vector<BlackMesh::BlackMesh<double>::BBX> cc_bbx;

	MeshSampling();
	~MeshSampling();

	//set mesh
	inline void setMesh(BlackMesh::BlackMesh<double> * msh_) {
		msh = msh_;
	}

	//sample  P=max{A/(diag*diag)*N, Pfix}
	void generate_sampling_points(const int total_points, const int minimal_points);

	inline void generate_area_and_other_things() {
		if (msh->GetNumComponents() == 0) {
			msh->mark_component_with_coherence();
		}
		assemble_data();
		assemble_bbx();
	}
};

