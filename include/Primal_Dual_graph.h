#pragma once
#include "BlackMesh.h"
#include "stdafx.h"

class Primal_Dual_graph {
public:
	struct Node {
		//the index of node should be the same as component id
		std::vector<int> incident_edges;//the edge id in this class;
										//int connected_component_id;
										//int id_in_compoennt;
	};

	struct Nonmanifold_edge {
		std::vector<int> incident_components;//should follow the order in the Edge.incident_face
		int original_edge_id;
		//int connected_component_id;
		//int degree;
	};

	std::vector<Node> node_list;
	std::vector<Nonmanifold_edge> edge_list;
	std::vector<std::vector<int>> non_manifold_curves;
	//std::vector<vector<int>> connected_component_list_edge;
	std::vector<vector<int>> connected_component_list_node;

	inline void clear() {
		node_list.clear();
		edge_list.clear();
		non_manifold_curves.clear();
		connected_component_list_node.clear();
		//mNumConnectedComponents = 0;
	}

	inline void build_connection_graph(BlackMesh::BlackMesh<double> *_msh)
	{
		if(_msh->GetNumComponents()==0)
		_msh->mark_component_with_coherence();

		node_list.clear();
		node_list.resize(_msh->GetNumComponents());

		_msh->update_edge_manifold_flag();
		for (int i = 0; i < _msh->GetNumEdges(); i++) {
			if (_msh->GetEdges()[i].is_nonmanifold_edge == 1) {
				Primal_Dual_graph::Nonmanifold_edge ed;
				for (int j = 0; j < _msh->GetEdges()[i].incident_faces.size(); j++) {
					int id =
						_msh->GetTriangles()[_msh->GetEdges()[i].incident_faces[j]].component_id;
					ed.incident_components.push_back(id);
					node_list[id].incident_edges.push_back(edge_list.size());
				}
				ed.original_edge_id = i;
				edge_list.push_back(ed);
			}
		}

	}


};

