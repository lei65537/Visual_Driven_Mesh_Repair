#pragma once
//#include "Public_header.h"

#include "Utils.h"
#include "BlackMesh.h"
#include "Primal_Dual_graph.h"
#include <set>
#include<Eigen/SparseCholesky>
#include <Eigen/Dense>
#include <igl/cotmatrix.h>

#define _ESP_FOR_FACE_CENTER_COMPARE 0.01

class eigen_vec_less {
public:
	//double eps;
	//eigen_vec_less(double arg_ = 1e-7) :eps(arg_) {}
	eigen_vec_less() {}
	bool operator()(const Eigen::Vector3d& v1, const Eigen::Vector3d & v2)const {
		return (v1[0] < v2[0]) || ((v1[0] == v2[0]) && (v1[1] < v2[1])) || ((v1[0] == v2[0]) && (v1[1] == v2[1]) && (v1[2] < v2[2]));
	}


};

template <class Real>
class Preprocess
{
private:
	//GEO::Mesh *_msh;
	double _merge_eps;
	//GEO::Mesh dual_mesh;
	double _shift_face_center;

	BlackMesh::BlackMesh<Real> *_msh;
	//BlackMesh::Dual_graph dual_mesh;

public:
	Preprocess( double merge_eps = 1e-10, double shift_face_center = 1e-2) :
		_msh(NULL), _merge_eps(merge_eps), _shift_face_center(shift_face_center) {
	}
	Preprocess(BlackMesh::BlackMesh<Real>& msh, double merge_eps=1e-7, double shift_face_center=1e-2) :
		_msh(&msh), _merge_eps(merge_eps), _shift_face_center(shift_face_center){
	}	
	Preprocess(BlackMesh::BlackMesh<Real>* msh, double merge_eps = 1e-7, double shift_face_center = 1e-2) :
		_msh(msh), _merge_eps(merge_eps), _shift_face_center(shift_face_center) {
	}
	~Preprocess() {}

	inline void set_mesh(BlackMesh::BlackMesh<Real> *msh) { _msh = msh; }

	//merge vertex
	void opt_vert()//;

	{
		int siz_v = _msh->GetNumVertices();
		int siz_f = _msh->GetNumTriangles();
		//merge identical vtx, use L2 norm here
		auto if_two_vec_identical = [=](Eigen::Vector3d &v1, Eigen::Vector3d &v2) {
			return (v1 - v2).norm() < _merge_eps;
		};

#ifdef _OUTPUT_INFO
		std::cout << "****Start Fetching Original Mesh Structure\n";
#endif
		std::vector<std::pair<Eigen::Vector3d, int>> vt_list_with_idx;
		std::vector<std::vector<int>> face_list;
		//get original mesh struct
		for (int i = 0; i < siz_v; i++) {
			Eigen::Vector3d tmp;
			//std::vector<double> pt = _msh->GetVertices()[i].pos;
			std::vector<double> pt;
			auto pp = _msh->GetVertices()[i].pos;
			Utils::to_double<Real>(pp,pt);
			tmp[0] = pt[0];
			tmp[1] = pt[1];
			tmp[2] = pt[2];
			vt_list_with_idx.push_back(std::make_pair(tmp, i));
		}
		for (int i = 0; i < siz_f; i++) {
			std::vector<int> tmp;
			for (int j = 0; j < _msh->GetTriangles()[i].vertices.size(); j++) {
				tmp.push_back(_msh->GetTriangles()[i].vertices[j]);
			}
			face_list.push_back(tmp);
		}




		std::vector<int> old_id_to_new_id(siz_v);
		old_id_to_new_id[vt_list_with_idx[0].second] = 0;

		std::vector<Eigen::Vector3d> new_vtx_list;
		new_vtx_list.push_back(vt_list_with_idx[0].first);


		for (int i = 0; i < vt_list_with_idx.size(); i++) {
			int currentvid = i;
			int push_flag = true;
			int hitid = -1;

			for (int j = 0; j < new_vtx_list.size(); j++) {
				if (if_two_vec_identical(new_vtx_list[j], vt_list_with_idx[i].first)) {
					push_flag = false;
					hitid = j;
					break;
				}
			}

			if (push_flag) {
				int siz = new_vtx_list.size();
				new_vtx_list.push_back(vt_list_with_idx[i].first);
				old_id_to_new_id[i] = siz;
			}
			else {
				old_id_to_new_id[i] = hitid;
			}
		}



#ifdef _OUTPUT_INFO
		std::cout << "****Building New Mesh \n";
#endif
		//reformulate vtx and face to GEO
		_msh->clear();
		for (int i = 0; i < new_vtx_list.size(); i++) {

			_msh->insert_vtx(std::vector<Real>{new_vtx_list[i][0],
				new_vtx_list[i][1], new_vtx_list[i][2]});
		}

		auto if_degenerate_face = [](std::vector<int> &tmp) {

			for (int i = 0; i < tmp.size(); i++) {
				for (int j = i + 1; j < tmp.size(); j++) {
					if (tmp[i] == tmp[j])
					{
						std::cout << "Hit";
						return true;
					}
				}
			}
			return false;
		};

	
		for (int i = 0; i < face_list.size(); i++) {
			
			std::vector<int> tmp = { old_id_to_new_id[face_list[i][0]],
				old_id_to_new_id[face_list[i][1]], old_id_to_new_id[face_list[i][2]] };
			if (!if_degenerate_face(tmp))
				_msh->insert_face(tmp);
		}
		_msh->update_mesh_properties();
	}

	void opt_vert_ext_on_pro(BlackMesh::BlackMesh<double> *prox_mesh, Real eps, bool zero_eps_flag = false, bool ext_to_double_norm = false) {

		int siz_v = prox_mesh->GetNumVertices();
		int siz_f = prox_mesh->GetNumTriangles();
		//merge identical vtx, use L2 norm here
		auto if_two_vec_identical = [=](Eigen::Vector3d &v1, Eigen::Vector3d &v2) {
			return (v1 - v2).norm() < _merge_eps;
		};

#ifdef _OUTPUT_INFO
		std::cout << "****Start Fetching Original Mesh Structure\n";
#endif
		std::vector<std::pair<Eigen::Vector3d, int>> vt_list_with_idx;
		std::vector<std::vector<int>> face_list;
		//get original mesh struct
		for (int i = 0; i < siz_v; i++) {
			Eigen::Vector3d tmp;
			std::vector<double> pt;
			auto pp = prox_mesh->GetVertices()[i].pos;
			tmp[0] = pp[0];
			tmp[1] = pp[1];
			tmp[2] = pp[2];
			vt_list_with_idx.push_back(std::make_pair(tmp, i));
		}
		for (int i = 0; i < siz_f; i++) {
			std::vector<int> tmp;
			for (int j = 0; j < prox_mesh->GetTriangles()[i].vertices.size(); j++) {
				tmp.push_back(prox_mesh->GetTriangles()[i].vertices[j]);
			}
			face_list.push_back(tmp);
		}




		std::vector<int> old_id_to_new_id(siz_v);
		old_id_to_new_id[vt_list_with_idx[0].second] = 0;

		std::vector<Eigen::Vector3d> new_vtx_list;
		new_vtx_list.push_back(vt_list_with_idx[0].first);


		for (int i = 0; i < vt_list_with_idx.size(); i++) {
			int currentvid = i;
			int push_flag = true;
			int hitid = -1;

			for (int j = 0; j < new_vtx_list.size(); j++) {
				if (if_two_vec_identical(new_vtx_list[j], vt_list_with_idx[i].first)) {
					push_flag = false;
					hitid = j;
					break;
				}
			}

			if (push_flag) {
				int siz = new_vtx_list.size();
				new_vtx_list.push_back(vt_list_with_idx[i].first);
				old_id_to_new_id[i] = siz;
			}
			else {
				old_id_to_new_id[i] = hitid;
			}
		}



#ifdef _OUTPUT_INFO
		std::cout << "****Building New Mesh \n";
#endif

		_msh->clear();

		for (int i = 0; i < new_vtx_list.size(); i++) {

			_msh->insert_vtx(std::vector<Real>{new_vtx_list[i][0],
				new_vtx_list[i][1], new_vtx_list[i][2]});
		}

		auto if_degenerate_face = [](std::vector<int> &tmp) {

			for (int i = 0; i < tmp.size(); i++) {
				for (int j = i + 1; j < tmp.size(); j++) {
					if (tmp[i] == tmp[j])
					{
						//std::cout << "Hit";
						return true;
					}
				}
			}
			return false;
		};


		for (int i = 0; i < face_list.size(); i++) {
			std::vector<int> tmp = { old_id_to_new_id[face_list[i][0]],
				old_id_to_new_id[face_list[i][1]], old_id_to_new_id[face_list[i][2]] };
			if (!if_degenerate_face(tmp))
				_msh->insert_face(tmp);
		}


		_msh->update_mesh_properties();
	}
	void opt_vert_ext(Real eps, bool zero_eps_flag=false, bool ext_to_double_norm=false)//;

	{
		typedef Eigen::Matrix<
			Real,
			3,
			1,
			Eigen::MatrixXd::IsRowMajor> Vec3ES;

		int siz_v = _msh->GetNumVertices();
		int siz_f = _msh->GetNumTriangles();
		//merge identical vtx, use L2 norm here
		auto ext_nm = [&](Vec3ES &v1, Vec3ES &v2) {
			Real x = v1[0] - v2[0];
			Real y = v1[1] - v2[1];
			Real z = v1[2] - v2[2];
			Real dis = x * x
				+ y * y 
				+ z  * z ;
			return dis;
		};
		auto if_two_vec_identical = [=](Vec3ES &v1, Vec3ES &v2, const Real &eps2, const bool &zero_eps_flag,  const bool &double_flag
			/*, bool outputFlag=false*/) {
			if (!double_flag) {
				if (!zero_eps_flag) {
					Real dis = ext_nm(v1, v2);
					return dis <= eps2;
				}
				else {
					return (v1[0] == v2[0]) && (v1[1] == v2[1]) && (v1[2] == v2[2]);
				}
			}
			else {
				Eigen::Vector3d v1d, v2d;
				v1d << Utils::to_double(v1[0]),Utils::to_double(v1[1]),Utils::to_double(v1[2]);
				v2d << Utils::to_double(v2[0]),Utils::to_double(v2[1]),Utils::to_double(v2[2]);
				return (v1d - v2d).squaredNorm() <= Utils::to_double(eps2);
			}
		};
		
#ifdef _OUTPUT_INFO
		std::cout << "****Start Fetching Original Mesh Structure\n";
#endif
		std::vector<std::pair<Vec3ES, int>> vt_list_with_idx;
		std::vector<std::vector<int>> face_list;
		//get original mesh struct
		for (int i = 0; i < siz_v; i++) {
			Vec3ES tmp;
			//std::vector<double> pt = _msh->GetVertices()[i].pos;
			std::vector<double> pt;
			auto pp = _msh->GetVertices()[i].pos;
			tmp[0] = pp[0];
			tmp[1] = pp[1];
			tmp[2] = pp[2];
			vt_list_with_idx.push_back(std::make_pair(tmp, i));
		}
		for (int i = 0; i < siz_f; i++) {
			std::vector<int> tmp;
			for (int j = 0; j < _msh->GetTriangles()[i].vertices.size(); j++) {
				tmp.push_back(_msh->GetTriangles()[i].vertices[j]);
			}
			face_list.push_back(tmp);
		}




		std::vector<int> old_id_to_new_id(siz_v);
		old_id_to_new_id[vt_list_with_idx[0].second] = 0;

		std::vector<Vec3ES> new_vtx_list;
		new_vtx_list.push_back(vt_list_with_idx[0].first);

		eps = eps*eps;
		for (int i = 0; i < vt_list_with_idx.size(); i++) {
			int currentvid = i;
			int push_flag = true;
			int hitid = -1;

			for (int j = 0; j < new_vtx_list.size(); j++) {

				if (if_two_vec_identical(new_vtx_list[j], vt_list_with_idx[i].first, eps, zero_eps_flag, ext_to_double_norm)) {
					push_flag = false;
					hitid = j;
					break;
				}
			}

			if (push_flag) {
				int siz = new_vtx_list.size();
				new_vtx_list.push_back(vt_list_with_idx[i].first);
				old_id_to_new_id[i] = siz;
			}
			else {
				old_id_to_new_id[i] = hitid;
			}
		}



#ifdef _OUTPUT_INFO
		std::cout << "****Building New Mesh \n";
#endif

		_msh->clear();

		for (int i = 0; i < new_vtx_list.size(); i++) {
			_msh->insert_vtx(std::vector<Real>{new_vtx_list[i][0],
				new_vtx_list[i][1], new_vtx_list[i][2]});
		}

		auto if_degenerate_face = [](std::vector<int> &tmp) {

			for (int i = 0; i < tmp.size(); i++) {
				for (int j = i + 1; j < tmp.size(); j++) {
					if (tmp[i] == tmp[j])
					{
						std::cout << "Hit";
						return true;
					}
				}
			}
			return false;
		};


		for (int i = 0; i < face_list.size(); i++) {
			std::vector<int> tmp = { old_id_to_new_id[face_list[i][0]],
				old_id_to_new_id[face_list[i][1]], old_id_to_new_id[face_list[i][2]] };
			if (!if_degenerate_face(tmp))
				_msh->insert_face(tmp);
		}

		_msh->update_mesh_properties();
	}


	//generate dual graph
	void generate_dual_graph(BlackMesh::Dual_graph<Real> &dual_mesh)//;

	{
		dual_mesh.clear();

		//face id = node id
		for (int i = 0; i < _msh->GetNumTriangles(); i++) {
			std::vector<double> ct = { 0, 0, 0 };
			for (int j = 0; j < _msh->GetTriangles()[i].vertices.size(); j++) {
				std::vector<double> pt = _msh->GetVertices()[_msh->GetTriangles()[i].vertices[j]].pos;
				ct[0] += pt[0];
				ct[1] += pt[1];
				ct[2] += pt[2];
			}
			ct[0] = ct[0] / (double)_msh->GetTriangles()[i].vertices.size();
			ct[1] = ct[1] / (double)_msh->GetTriangles()[i].vertices.size();
			ct[2] = ct[2] / (double)_msh->GetTriangles()[i].vertices.size();

			std::vector<double> n_vt = { 0,0,0 };

			Eigen::Vector3d displacement = Eigen::Vector3d::Random();
			displacement.normalize();

			n_vt[0] = ct[0] + displacement[0] * _shift_face_center;
			n_vt[1] = ct[1] + displacement[1] * _shift_face_center;
			n_vt[2] = ct[2] + displacement[2] * _shift_face_center;

			dual_mesh.insert_vtx(n_vt);
		}
		std::cout << "Dual Mesh has " << dual_mesh.GetNumVertices() << " vertices\n";
		//go through all faces, make the adjacent relationship beween faces to edge of the graph, for face i only count face j when j > i
		std::vector<std::pair<int, int>> debug_edge;
		for (int i = 0; i < _msh->GetNumTriangles(); i++) {
			for (int j = 0; j < _msh->GetTriangles()[i].incident_edges.size(); j++) {
				int inc_edge = _msh->GetTriangles()[i].incident_edges[j];
				for (int k = 0; k < _msh->GetEdges()[inc_edge].incident_faces.size(); k++) {
					//for each edge go through all incident faces
					if (_msh->GetEdges()[inc_edge].incident_faces[k] != i) {
						//not idential face, base face i, incident face _msh->GetEdges()[inc_edge].incident_faces[k], make sure j > i
						dual_mesh.insert_edge(i, _msh->GetEdges()[inc_edge].incident_faces[k], inc_edge);
					}
				}
			}
		}


	}

	//delete duplicate faces, if two faces share all three edges, delete one
	void delete_duplicate_faces()//;

	{

		double ave_edge_lenth = _msh->GetAverageEdge();
		assert(ave_edge_lenth != 0);

		auto if_two_face_share_same_three_edges = [](BlackMesh::BlackMesh<Real>::Triangle &fc1, BlackMesh::BlackMesh<Real>::Triangle &fc2) {
			int hit_id = 0;
			for (int i = 0; i < fc1.incident_edges.size(); i++) {
				for (int j = 0; j < fc2.incident_edges.size(); j++) {
					if (fc1.incident_edges[i] == fc2.incident_edges[j]) {
						hit_id++;
						break;
					}
				}
			}

			if (hit_id == 3) return true;
			else return false;
		};

		std::map<int, int> check_list;//faceid,0-keep 1-del

		for (int i = 0; i < _msh->GetNumEdges(); i++) {
			if (_msh->GetEdges()[i].incident_faces.size() > 1)//check if edge is non manifold// boundary edge 
			{

				BlackMesh::BlackMesh<Real>::Edge ed = _msh->GetEdges()[i];
				for (int j = 0; j < ed.incident_faces.size(); j++) {
					int base_fc_id = ed.incident_faces[j];
					BlackMesh::BlackMesh<Real>::Triangle base_fc = _msh->GetTriangles()[base_fc_id];
					Eigen::Vector3d base_fc_ct(base_fc.face_center.data());


					for (int k = j + 1; k < ed.incident_faces.size(); k++) {
						int compare_fc_id = ed.incident_faces[k];
						//if (check_list[compare_fc_id] == true) continue;

						BlackMesh::BlackMesh<Real>::Triangle compare_fc = _msh->GetTriangles()[compare_fc_id];
						Eigen::Vector3d compare_fc_ct(compare_fc.face_center.data());

						//tell if two face have close center
						if ((base_fc_ct - compare_fc_ct).norm() < _ESP_FOR_FACE_CENTER_COMPARE*ave_edge_lenth) {
							//double check if the three edges are same
							bool flg = if_two_face_share_same_three_edges(compare_fc, base_fc);


							if (flg) {
								auto itbase = check_list.find(base_fc_id);
								auto itkeep = check_list.find(compare_fc_id);

								if (itbase == check_list.end() && itkeep == check_list.end()) {
									//both not visited

									check_list.insert(std::make_pair(base_fc_id, 0));
									check_list.insert(std::make_pair(compare_fc_id, 1));

								}
								else if (itbase == check_list.end() && itkeep != check_list.end()) {
									check_list.insert(std::make_pair(base_fc_id, 1));
								}
								else if (itbase != check_list.end() && itkeep == check_list.end()) {
									check_list.insert(std::make_pair(compare_fc_id, 1));
								}

							}

						}

					}
				}


			}
		}


		//construct a new mesh
		std::vector<BlackMesh::BlackMesh<Real>::Triangle> tri_copy = _msh->GetTrianglesCopy();

		_msh->clear_but_keep_vtx();

		for (int i = 0; i < tri_copy.size(); i++) {
			auto it = check_list.find(i);
			if (it == check_list.end() || (it->second) == 0) {

				_msh->insert_face(tri_copy[i].vertices);
			}
		}

		_msh->update_mesh_properties();

	}

	//
	void get_max_connected_mesh_list(std::vector<BlackMesh::BlackMesh<Real>*> &output_msh_list)//;

	{
		Primal_Dual_graph gh;
		gh.build_connection_graph(_msh);

		//flood from node
		std::vector<bool> visit_flag;
		visit_flag.resize(gh.node_list.size(), false);
		std::vector<std::set<int>> grouped_node;

		for (int m = 0; m < visit_flag.size(); m++) {
			if (visit_flag[m])continue;

			std::queue<int> list;
			list.push(m);
			grouped_node.push_back(std::set<int>{m});
			visit_flag[0] = true;

			int gpid = grouped_node.size();
			gpid--;

			while (list.size() != 0) {
				int current_node = list.front();
				list.pop();

				auto nh = gh.node_list[current_node];
				for (int i = 0; i < nh.incident_edges.size(); i++) {
					int ed = nh.incident_edges[i];
					for (int j = 0; j < gh.edge_list[ed].incident_components.size(); j++) {
						int inciid = gh.edge_list[ed].incident_components[j];
						if (visit_flag[inciid] == false) {
							visit_flag[inciid] = true;
							list.push(inciid);
							grouped_node[gpid].insert(inciid);
						}
					}
				}

			}
		}

		//split the mesh;
		output_msh_list.clear();
		for (int mm = 0; mm < grouped_node.size(); mm++) {
			//output_msh_list.push_back(BlackMesh::BlackMesh());
			auto ptr = new BlackMesh::BlackMesh<Real>();
			ptr->copy_vertice_from(_msh);

			for (int i = 0; i < _msh->GetNumTriangles(); i++) {
				auto it = grouped_node[mm].find(_msh->GetTriangles()[i].component_id);
				if (it != grouped_node[mm].end()) {
					auto vtx = (_msh->GetTriangles()[i].vertices);
					ptr->insert_face(vtx);
				}
			}
			output_msh_list.push_back(ptr);
		}
	}

	//check if all faces are not single
	bool check_duplicate_face()//;

	{
		for (int i = 0; i < _msh->GetNumTriangles(); i++) {
			auto vtxbase = _msh->GetTriangles()[i].vertices;

			std::set<int> setvtxbase;
			setvtxbase.insert(vtxbase.begin(), vtxbase.end());

			for (int j = i + 1; j < _msh->GetNumTriangles(); j++) {

				auto vtxquery = _msh->GetTriangles()[j].vertices;
				if (setvtxbase.find(vtxquery[0]) != setvtxbase.end() &&
					setvtxbase.find(vtxquery[1]) != setvtxbase.end() &&
					setvtxbase.find(vtxquery[2]) != setvtxbase.end())
				{
					return false;
				}
			}
		}

		return true;
	}

	//make the cc coherent as much as possible, only on 2-edge
	void make_cc_coherent(std::vector<double> &posScore, std::vector<double> &negScore)//;

	{
		Primal_Dual_graph gh;

		BlackMesh::BlackMesh<double> shadowMesh;
		shadowMesh.construct_from<Real>(*_msh);
		gh.build_connection_graph(&shadowMesh);


		std::vector<bool> visited_flag, flip_flag;
		visited_flag.resize(gh.node_list.size(), false);
		flip_flag.resize(gh.node_list.size(), false);

		std::vector<std::pair<int, int>> cc_nfc;//number of faces in a cc
		cc_nfc.reserve(gh.node_list.size());
		for (int i = 0; i < gh.node_list.size(); i++) {
			cc_nfc.push_back(std::make_pair(i, 0));
		}
		for (int i = 0; i < _msh->GetNumTriangles(); i++) {
			cc_nfc[_msh->GetTriangles()[i].component_id].second += 1;
		}

		std::sort(cc_nfc.begin(), cc_nfc.end(), [](std::pair<int, int>& a, std::pair<int, int>& b) {
			return a.second < b.second;
		});

		while (cc_nfc.size() != 0) {
			int siz = cc_nfc.size();
			int i = cc_nfc[siz - 1].first;
			cc_nfc.pop_back();

			if (visited_flag[i]) continue;

			std::queue<int> list;
			list.push(i);
			visited_flag[i] = true;

			if (negScore[i] > posScore[i]) {
				flip_flag[i] = true;
			}

			while (list.size() != 0) {
				int current_cc = list.front();
				list.pop();

				for (int j = 0; j < gh.node_list[current_cc].incident_edges.size(); j++) {
					int edgeid = gh.node_list[current_cc].incident_edges[j];
					auto eh = gh.edge_list[edgeid];
					if (eh.incident_components.size() != 2)continue;

					//cci -> ccj through a 2-edge
					int ccj = eh.incident_components[0] == current_cc ?
						eh.incident_components[1] :
						eh.incident_components[0];

					if (visited_flag[ccj] == false) {
						bool dir_flg = !flip_flag[current_cc];//dir_flg is the direction ccj if connect; true means flip

						bool go_on_flag = true;
						if (dir_flg) {
							//ccj needs to flip to connect
							//the score of flip must larger than keep its original dir
							if (negScore[ccj] < posScore[ccj]) {
								go_on_flag = false;
							}
						}
						if (go_on_flag) {
							visited_flag[ccj] = true;
							flip_flag[ccj] = !flip_flag[current_cc];
							list.push(ccj);
						}
					}
				}
			}
		}

		//updated mesh
		for (int i = 0; i < _msh->GetNumTriangles(); i++) {
			if (flip_flag[_msh->GetTriangles()[i].component_id]) {
				auto vtxcp = _msh->GetTriangles()[i].vertices;
				_msh->GetTrianglesEditable()[i].vertices[1] = vtxcp[2];
				_msh->GetTrianglesEditable()[i].vertices[2] = vtxcp[1];
			}
		}

		_msh->update_mesh_properties();

	}

	void make_cc_coherent(std::vector<double> &posScore, std::vector<double> &negScore, 
		BlackMesh::BlackMesh<double> &shadowMesh)//;

	{
		Primal_Dual_graph gh;


		gh.build_connection_graph(&shadowMesh);

		//gh.stack_non_manifold_curves(_msh);

		std::vector<bool> visited_flag, flip_flag;
		visited_flag.resize(gh.node_list.size(), false);
		flip_flag.resize(gh.node_list.size(), false);

		std::vector<std::pair<int, int>> cc_nfc;//number of faces in a cc
		cc_nfc.reserve(gh.node_list.size());
		for (int i = 0; i < gh.node_list.size(); i++) {
			cc_nfc.push_back(std::make_pair(i, 0));
		}
		for (int i = 0; i < _msh->GetNumTriangles(); i++) {
			cc_nfc[_msh->GetTriangles()[i].component_id].second += 1;
		}

		std::sort(cc_nfc.begin(), cc_nfc.end(), [](std::pair<int, int>& a, std::pair<int, int>& b) {
			return a.second < b.second;
		});

		while (cc_nfc.size() != 0) {
			int siz = cc_nfc.size();
			int i = cc_nfc[siz - 1].first;
			cc_nfc.pop_back();

			if (visited_flag[i]) continue;

			std::queue<int> list;
			list.push(i);
			visited_flag[i] = true;

			if (negScore[i] > posScore[i]) {
				flip_flag[i] = true;
			}

			while (list.size() != 0) {
				int current_cc = list.front();
				list.pop();

				for (int j = 0; j < gh.node_list[current_cc].incident_edges.size(); j++) {
					int edgeid = gh.node_list[current_cc].incident_edges[j];
					auto eh = gh.edge_list[edgeid];
					if (eh.incident_components.size() != 2)continue;

					//cci -> ccj through a 2-edge
					int ccj = eh.incident_components[0] == current_cc ?
						eh.incident_components[1] :
						eh.incident_components[0];

					if (visited_flag[ccj] == false) {
						bool dir_flg = !flip_flag[current_cc];//dir_flg is the direction ccj if connect; true means flip

						bool go_on_flag = true;
						if (dir_flg) {
							//ccj needs to flip to connect
							//the score of flip must larger than keep its original dir
							if (negScore[ccj] < posScore[ccj]) {
								go_on_flag = false;
							}
						}
						if (go_on_flag) {
							visited_flag[ccj] = true;
							flip_flag[ccj] = !flip_flag[current_cc];
							list.push(ccj);
						}
					}
				}
			}
		}

		//updated mesh
		for (int i = 0; i < _msh->GetNumTriangles(); i++) {
			if (flip_flag[_msh->GetTriangles()[i].component_id]) {
				auto vtxcp = _msh->GetTriangles()[i].vertices;
				_msh->GetTrianglesEditable()[i].vertices[1] = vtxcp[2];
				_msh->GetTrianglesEditable()[i].vertices[2] = vtxcp[1];

				auto vtxcps = shadowMesh.GetTriangles()[i].vertices;
				shadowMesh.GetTrianglesEditable()[i].vertices[1] = vtxcps[2];
				shadowMesh.GetTrianglesEditable()[i].vertices[2] = vtxcps[1];
			}
		}

		_msh->update_mesh_properties();
		shadowMesh.update_mesh_properties();
	}

	//shrink the boundary a little
	void shrink_boundary(double move_dis,double max_move_dis,  std::vector<std::pair<int,std::vector<double>>> *target_container=NULL)//;

	{
		/***********Useful Functions************/
		auto cross_prod = [](double *v1, double* v2, double* output) {
			output[0] = v1[1] * v2[2] - v1[2] * v2[1];
			output[1] = v1[2] * v2[0] - v1[0] * v2[2];
			output[2] = v1[0] * v2[1] - v1[1] * v2[0];
		};


		auto normalize = [](double * v) {
			double Length = sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
			if (Length < 1e-10) return 0;
			double invLength = (double)1 / Length;
			v[0] *= invLength;
			v[1] *= invLength;
			v[2] *= invLength;
			return 1;
		};


		auto modify_vtx_pos = [&](int eid1, int eid2) {
			int share_vid = -1;
			auto eh1 = _msh->GetEdges()[eid1];
			auto eh2 = _msh->GetEdges()[eid2];

			if (eh1.vertices.first == eh2.vertices.first) {
				share_vid = eh1.vertices.first;
			}
			else if (eh1.vertices.first == eh2.vertices.second) {
				share_vid = eh1.vertices.first;
			}
			else if (eh1.vertices.second == eh2.vertices.first) {
				share_vid = eh1.vertices.second;
			}
			else if (eh1.vertices.second == eh2.vertices.second) {
				share_vid = eh1.vertices.second;
			}
			else {
				//cannot be
			}

			auto tri1 = eh1.incident_faces[0];
			auto tri2 = eh2.incident_faces[0];

			double hei1[3];
			double hei2[3];

			double flag1 = _msh->compute_fake_height(tri1, eid1, hei1, 0, _SEARCH_NEIB_ON_DEGENERATE);
			double flag2 = _msh->compute_fake_height(tri2, eid2, hei2, 0, _SEARCH_NEIB_ON_DEGENERATE);

			if (flag1<0 || flag2<0) {
				std::cout << "Error";
#if defined(_RECORD_ERROR_LOGS)
				(*_log_off) << _file_under_processing << " ERROR Degenerate exist in boundary shrink\n";
				(*_log_off).flush();
#endif
			}

			double scale_Factor = sqrt(0.5*(flag1 + flag2));
			double move_Factor = scale_Factor*move_dis;
			if (move_Factor > max_move_dis)move_Factor = max_move_dis;

			if (target_container == NULL) {
				_msh->GetVerticesEditable()[share_vid].pos[0] += move_Factor*(hei1[0] + hei2[0]);
				_msh->GetVerticesEditable()[share_vid].pos[1] += move_Factor*(hei1[1] + hei2[1]);
				_msh->GetVerticesEditable()[share_vid].pos[2] += move_Factor*(hei1[2] + hei2[2]);
			}
			else {
				target_container->push_back(std::make_pair(share_vid,
					std::vector<double>{_msh->GetVerticesEditable()[share_vid].pos[0] + move_Factor*(hei1[0] + hei2[0]),
					_msh->GetVerticesEditable()[share_vid].pos[1] + move_Factor*(hei1[1] + hei2[1]),
					_msh->GetVerticesEditable()[share_vid].pos[2] + move_Factor*(hei1[2] + hei2[2])}));
			}
		};



		for (int i = 0; i < _msh->GetNumVertices(); i++) {
			auto vh = _msh->GetVertices()[i];
			std::vector<int> bound_edges;
			for (int j = 0; j < vh.incident_edges.size(); j++) {
				int edgeid = vh.incident_edges[j];
				if (_msh->GetEdges()[edgeid].incident_faces.size() == 1) {
					bound_edges.push_back(edgeid);
				}
			}

			if (bound_edges.size() == 2) {
				modify_vtx_pos(bound_edges[0], bound_edges[1]);
			}
			else if (bound_edges.size() == 0) {

			}
			else {
				std::cout << "error";
#if defined(_RECORD_ERROR_LOGS)
				(*_log_off) << _file_under_processing << " ERROR boundary edge degree weird\n";
				(*_log_off).flush();
#endif
			}
		}

		_msh->update_mesh_properties();
	}

	void shrink_boundary_laplacian(double move_dis, double lapweight=0.001, double boundweight=0.1, double posweight=1)//;

	{
		
	}

	void shrink_boundary_curve_lap(double move_dis, double max_move_dis, double lapweight = 0.05, double boundweight = 0.1, double posweight = 1)//;

	{

		typedef Eigen::SparseMatrix<double> SpMat;
		typedef Eigen::Triplet<double> Trip;

		lapweight = sqrt(lapweight);
		boundweight = sqrt(boundweight);
		posweight = sqrt(posweight);
		auto if_is_bd_vtx = [&](int vid, std::vector<int> &neiid) {
			int counter = 0;
			for (int i = 0; i < _msh->GetVertices()[vid].incident_edges.size(); i++) {
				int eid = _msh->GetVertices()[vid].incident_edges[i];
				if (_msh->GetEdges()[eid].incident_faces.size() == 1) {
					int anid = _msh->GetEdges()[eid].vertices.first == vid ?
						_msh->GetEdges()[eid].vertices.second :
						_msh->GetEdges()[eid].vertices.first;
					neiid.push_back(anid);
					counter++;
				}
			}
			return counter == 2;
		};

		std::vector<Trip> trilist;
		int cons_counter = 0;
		std::vector<std::pair<int, std::vector<double>>> target_container;
		shrink_boundary(move_dis, max_move_dis, &target_container);

		std::vector<bool> visited_flag;
		visited_flag.resize(_msh->GetNumVertices(), false);

		Eigen::VectorXd RHSX(_msh->GetNumVertices());
		Eigen::VectorXd RHSY(_msh->GetNumVertices());
		Eigen::VectorXd RHSZ(_msh->GetNumVertices());

		for (int i = 0; i < target_container.size(); i++) {
			int vid = target_container[i].first;
			trilist.push_back(Trip(cons_counter, vid, boundweight));
			cons_counter++;
			//LtL.coeffRef(vid, vid) += boundweight*boundweight;
			RHSX(vid) = boundweight*boundweight*target_container[i].second[0];
			RHSY(vid) = boundweight*boundweight*target_container[i].second[1];
			RHSZ(vid) = boundweight*boundweight*target_container[i].second[2];
			visited_flag[vid] = true;
		}
		for (int i = 0; i < _msh->GetNumVertices(); i++) {
			if (visited_flag[i]) {
			}
			else {
				trilist.push_back(Trip(cons_counter, i, posweight));
				cons_counter++;
				//LtL.coeffRef(i, i) += posweight*posweight;
				RHSX(i) = posweight*posweight*_msh->GetVertices()[i].pos[0];
				RHSY(i) = posweight*posweight*_msh->GetVertices()[i].pos[1];
				RHSZ(i) = posweight*posweight*_msh->GetVertices()[i].pos[2];
			}
		}

		for (int i = 0; i < _msh->GetNumVertices(); i++) {
			std::vector<int> neiid;
			bool flg = if_is_bd_vtx(i, neiid);
			if (flg) {
				trilist.push_back(Trip(cons_counter, i, 1 * lapweight));
#ifdef _OUTPUT_DEBUG_LOGS
				if (neiid.size() != 2) {
					std::cout << "cannot be";
#if defined(_RECORD_ERROR_LOGS)
					(*_log_off) << _file_under_processing << " ERROR boundary vtx is not degree 2 bdedge\n";
					(*_log_off).flush();
#endif
				}
#endif

				for (int j = 0; j < neiid.size(); j++) {
					trilist.push_back(Trip(cons_counter, neiid[j], -0.5*lapweight));
				}
				cons_counter++;
			}
		}

		SpMat L(cons_counter, _msh->GetNumVertices()), LtL;
		L.setFromTriplets(trilist.begin(), trilist.end());
		LtL = L.transpose()*L;

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(LtL);
		Eigen::VectorXd Rx = solver.solve(RHSX);
		Eigen::VectorXd Ry = solver.solve(RHSY);
		Eigen::VectorXd Rz = solver.solve(RHSZ);

		for (int i = 0; i < _msh->GetNumVertices(); i++) {
			_msh->GetVerticesEditable()[i].pos[0] = Rx(i);
			_msh->GetVerticesEditable()[i].pos[1] = Ry(i);
			_msh->GetVerticesEditable()[i].pos[2] = Rz(i);
		}

		_msh->update_mesh_properties();
	}
};
