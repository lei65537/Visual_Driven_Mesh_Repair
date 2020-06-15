#include "stdafx.h"
#include "Optimizer.h"


Optimizer::Optimizer()
{

}

Optimizer::~Optimizer()
{

}

int Optimizer::if_face_coherant_based_on_edge(int vid_start, int vid_end, 
	const BlackMesh::BlackMesh<double>::Triangle &tr1, const BlackMesh::BlackMesh<double>::Triangle &tr2) {
	//1-coherant
	std::vector<int> vtx1 = tr1.vertices;
	std::vector<int> vtx2 = tr2.vertices;

	int combed_vid_s = -1;
	int combed_vid_e = -1;

	for (int i = 0; i < vtx1.size(); i++) {
		if (vtx1[i] == vid_start) {
			if (vtx1[(i + 1) % 3] == vid_end) {
				combed_vid_s = vid_start;
				combed_vid_e = vid_end;
			}
			else if (vtx1[(i + 2) % 3] == vid_end) {
				combed_vid_s = vid_end;
				combed_vid_e = vid_start;
			}
		}
	}

	for (int i = 0; i < vtx2.size(); i++) {
		if (vtx2[i] == combed_vid_e&&vtx2[(i + 1) % 3] == combed_vid_s) {
			return 1;
		}
		else if (vtx2[i] == combed_vid_s&&vtx2[(i + 1) % 3] == combed_vid_e) {
			return 0;
		}
	}
	return -1;
};

void Optimizer::turn_conf_to_stacked(std::vector<std::vector<double>> &xin, std::vector<int> &xout) {
	xout.clear();

	std::vector<int> stacked_curve_x;
	for (int i = 0; i < xin.size(); i++) {
		for (int j = 0; j < xin[i].size(); j++) {
			if (abs(xin[i][j] - 1.0) < 1e-7)
			{
				stacked_curve_x.push_back(j);
				break;
			}
		}
	}
	xout.resize(_gh.edge_list.size(), -1);
	for (int m = 0; m < _gh.non_manifold_curves.size(); m++) {
		for (int k = 0; k < _gh.non_manifold_curves[m].size(); k++) {
			int edgeid = _gh.non_manifold_curves[m][k];
			xout[edgeid] = stacked_curve_x[m];
		}
	}

}

//init
void Optimizer::generate_transformation_matrices()
{
	configuration_list_all_degrees.clear();
	dual_edge_list_all_degrees.clear();

	trans_conf2dual_eign.clear();

	int max_degree = _msh->GetMaxEdgeDegree();
	for (int i = 2; i <= max_degree; i++) {
		//from degree 2 edge to the max degree

		std::vector<int> idxrow_conf2dual, idxcol_conf2dual;

		//generate a directional dual edge list
		std::vector<dual_edge_st> dual_edge_list;
		for (int m = 0; m < i; m++) {
			for (int n = m + 1; n < i; n++) {
				dual_edge_list.push_back(dual_edge_st(m, n, 0));
				dual_edge_list.push_back(dual_edge_st(m, n, 1));
			}
		}

		int genuine_dual_siz = dual_edge_list.size();

		//push virtual dual edge, used for orientation of un-connected vtx
		for (int m = 0; m < i; m++) {
			dual_edge_list.push_back(dual_edge_st(m, -1, 0));
			dual_edge_list.push_back(dual_edge_st(m, -1, 1));
		}

		//list all configurations
		std::vector<std::vector<int>> confs_last_time, confs_all;
		for (int m = i; m > 0; m = m - 2) {
			//list the configuration with i-m+2 edge connected
			add_one_dual_edge_on_current_configurations(confs_last_time, dual_edge_list);
			confs_all.insert(confs_all.end(), confs_last_time.begin(), confs_last_time.end());
		}
		//and one all zero configuration should be added
		confs_all.push_back(std::vector<int>{});

		std::vector<std::vector<int>> output;
		split_unorientated_cases(confs_all, dual_edge_list, i, genuine_dual_siz, output);
		confs_all.clear();
		confs_all.insert(confs_all.begin(), output.begin(), output.end());

		int siz_conf2dual = (confs_all.size());
		//int siz_conf2face=(dual_edge_list.size());
		int siz_conf2face = genuine_dual_siz;

		configuration_list_all_degrees.push_back(confs_all);
		dual_edge_list_all_degrees.push_back(dual_edge_list);

		//assemble matrix
		std::vector<Tri> trilist;
		for (int m = 0; m < confs_all.size(); m++) {
			std::vector<int> conf = confs_all[m];
			for (int s = 0; s < conf.size(); s++) {
				if (conf[s] >= genuine_dual_siz) continue;
				idxrow_conf2dual.push_back(conf[s]);
				idxcol_conf2dual.push_back(m);
				trilist.push_back(Tri(conf[s], m, 1.0));
			}
		}

		SpMat tmpsp(siz_conf2face, siz_conf2dual);//#duals #confs
		tmpsp.setFromTriplets(trilist.begin(), trilist.end());
		trans_conf2dual_eign.push_back(tmpsp);


	}

}

int  Optimizer::get_dual_edge_id(int fi, int fj, int degree, int direction) {
	//i and j are the index in the adjacent list, not face id
	//make sure j > i
	//direction-0/1
	if (fi == 0) return fj;
	if (fi == degree - 1) return -1;//this should not happen

	int counter_s = degree - fi;
	int counter_e = degree - 1;
	int term_counter = fi;
	int residul = fj - fi;

	return ((counter_e + counter_s)*term_counter / 2 + residul) * 2 + direction;
}

void  Optimizer::add_one_dual_edge_on_current_configurations(std::vector<std::vector<int>> &confs_res, std::vector<dual_edge_st> &dual_edge_list) {
	std::vector<std::vector<int>> confs = confs_res;
	confs_res.clear();
	if (confs.size() == 0) {
		//initialization
		for (int i = 0; i < dual_edge_list.size(); i++) {
			if (dual_edge_list[i].fj != -1)
				confs_res.push_back(std::vector<int>{i});
		}
		return;
	}

	//go through all existing cases in confs and add one more
	for (int i = 0; i < confs.size(); i++) {
		std::vector<int> current_conf = confs[i];
		//fetch max dual edge id
		int max_dual_id = -1;
		for (int j = 0; j < current_conf.size(); j++) {
			if (max_dual_id < current_conf[j])max_dual_id = current_conf[j];
		}

		//add one more, start from maxid,
		for (int j = max_dual_id + 1; j < dual_edge_list.size(); j++) {
			//try add jth edge, this edge should not have same face id as existing ones
			if (dual_edge_list[j].fj == -1)continue;
			bool conflict_flag = false;
			for (int m = 0; m < current_conf.size(); m++) {
				int testid = current_conf[m];
				if (dual_edge_list[j].fi == dual_edge_list[testid].fi
					|| dual_edge_list[j].fi == dual_edge_list[testid].fj
					|| dual_edge_list[j].fj == dual_edge_list[testid].fi
					|| dual_edge_list[j].fj == dual_edge_list[testid].fj) {
					conflict_flag = true;
					break;
				}
			}

			if (conflict_flag == false) {
				//not conflict, should add
				std::vector<int> tmp = current_conf;
				tmp.push_back(j);
				confs_res.push_back(tmp);
			}
		}
	}
}



void Optimizer::split_unorientated_cases(std::vector<std::vector<int>>& confs_res, std::vector<dual_edge_st>& dual_edge_list, int degree, int dual_siz, std::vector<std::vector<int>>& output)
{
	//std::vector<std::vector<int>> output;

	for (int i = 0; i < confs_res.size(); i++) {
		std::vector<bool> face_connected_flag;
		face_connected_flag.resize(degree, false);

		std::vector<std::vector<int>> this_conf;
		this_conf.push_back(confs_res[i]);
		for (int j = 0; j < confs_res[i].size(); j++) {
			int fidi = dual_edge_list[confs_res[i][j]].fi;
			int fidj = dual_edge_list[confs_res[i][j]].fj;
			face_connected_flag[fidi] = true;
			face_connected_flag[fidj] = true;
		}

		for (int j = 0; j < face_connected_flag.size(); j++) {
			if (face_connected_flag[j]) continue;
			//double the conf and branch element j
			int siz = this_conf.size();
			this_conf.insert(this_conf.begin(), this_conf.begin(), this_conf.end());
			for (int m = 0; m < siz; m++) {
				this_conf[m].push_back(dual_siz + 2 * j + 0);
			}
			for (int m = siz; m < 2 * siz; m++) {
				this_conf[m].push_back(dual_siz + 2 * j + 1);
			}
		}

		output.insert(output.end(), this_conf.begin(), this_conf.end());
	}
}



void Optimizer::solve_ccdir_on_current_sol(std::vector<vector<double>>& sol, std::map<int, int> *node, std::vector<int>& ccdir)
{


	//assemble fixed unary tag
	std::map<int, bool> node_induced_dir;
	if (node != NULL) {
		for (auto it = node->begin(); it != node->end(); it++) {
			int picked_curveid = it->first;
			int picked_confid = it->second;

			int degree = _gh.edge_list
				[_gh.non_manifold_curves[picked_curveid][0]].incident_components.size();

			for (int i = 0;
				i < configuration_list_all_degrees[degree - 2][picked_confid].size();
				i++) {

				int dualedgeid = configuration_list_all_degrees[degree - 2][picked_confid][i];

				auto dualedge = dual_edge_list_all_degrees[degree - 2][dualedgeid];

				if (dualedge.fj == -1) {
					node_induced_dir.insert(std::make_pair(_gh.edge_list
						[_gh.non_manifold_curves[picked_curveid][0]].incident_components[dualedge.fi], dualedge.dir == 0));
					continue;
				}

				auto eh = _gh.edge_list
					[_gh.non_manifold_curves[picked_curveid][0]];

				int ccidi = eh.incident_components[dualedge.fi];
				int ccidj = eh.incident_components[dualedge.fj];

				if (ccidi == ccidj) {
					node_induced_dir.insert(std::make_pair(ccidi, dualedge.dir == 0));

				}
				else {
					auto pit = binary_score.find(std::make_pair(picked_curveid, dualedgeid));
					node_induced_dir.insert(std::make_pair(eh.incident_components[pit->second.inidi], pit->second.diri == 0));
					node_induced_dir.insert(std::make_pair(eh.incident_components[pit->second.inidj], pit->second.dirj == 0));
				}


			}
		}
	}

	//accumulate unary and binary score
	std::vector<double> unary_term;
	unary_term.resize(2 * _gh.node_list.size(), 0.0);
	std::map<std::pair<int, int>, std::vector<double>> binary_term;//ccidi,ccidj, 00,01,10,11;

	for (int i = 0; i < unary_score.size(); i++) {
		//accumulate prob of conf
		{
			double xprob = 0.0;
			for (auto it = unary_score[i].confid_pos.begin();
				it != unary_score[i].confid_pos.end(); it++) {

				xprob += sol[it->first][it->second];
			}

			unary_term[2 * i] = unary_score[i].pos_score*xprob*unary_score[i].weight;

		}

		{
			double xprob = 0.0;
			for (auto it = unary_score[i].confid_neg.begin();
				it != unary_score[i].confid_neg.end(); it++) {

				xprob += sol[it->first][it->second];
			}

			unary_term[2 * i + 1] = unary_score[i].neg_score*xprob*unary_score[i].weight;

		}
	}


	//unary from binary
	for (auto it = unary_score_from_binary.begin();
		it != unary_score_from_binary.end(); it++) {
		int curveid = it->first.first;
		int dualedgeid = it->first.second;

		int someedgeid = _gh.non_manifold_curves[curveid][0];
		int ccidi = _gh.edge_list[someedgeid].incident_components[it->second.inidi];
		int ccidj = _gh.edge_list[someedgeid].incident_components[it->second.inidj];

#ifdef _OUTPUT_DEBUG_LOGS
		if (ccidi != ccidj) {
			std::cout << "Error\n";
#if defined(_RECORD_ERROR_LOGS)
			(*_log_off) << _file_under_processing << " ERROR 6\n";
			(*_log_off).flush();
#endif
		}
#endif
		double xprob = 0;
		for (auto ti = it->second.confid.begin();
			ti != it->second.confid.end(); ti++) {
			int confid = *ti;
			xprob += sol[curveid][confid];
		}


		double sc = it->second.weight*it->second.score*xprob;


		unary_term[2 * ccidi + it->second.diri] += sc;
	}

	//binary term
	for (auto it = binary_score.begin();
		it != binary_score.end(); it++) {

		int curveid = it->first.first;
		int dualedgeid = it->first.second;

		int someedgeid = _gh.non_manifold_curves[curveid][0];
		int ccidi = _gh.edge_list[someedgeid].incident_components[it->second.inidi];
		int ccidj = _gh.edge_list[someedgeid].incident_components[it->second.inidj];

#ifdef _OUTPUT_DEBUG_LOGS
		if (ccidi == ccidj) {
			std::cout << "Error\n";
#if defined(_RECORD_ERROR_LOGS)
			(*_log_off) << _file_under_processing << " ERROR 7\n";
			(*_log_off).flush();
#endif
		}
#endif

		double xprob = 0;
		for (auto ti = it->second.confid.begin();
			ti != it->second.confid.end(); ti++) {
			int confid = *ti;
			xprob += sol[curveid][confid];
		}


		double sc = it->second.weight*it->second.score*xprob;

		int diri = it->second.diri;
		int dirj = it->second.dirj;
		if (ccidi < ccidj) {
			int idx = diri * 2 + dirj;
			auto cit = binary_term.insert(std::make_pair(std::make_pair(ccidi, ccidj), std::vector<double>(4)));
			if (cit.second == true) {
				cit.first->second[0] = 0.0;
				cit.first->second[1] = 0.0;
				cit.first->second[2] = 0.0;
				cit.first->second[3] = 0.0;
			}
			cit.first->second[idx] += sc;
		}
		else {
			int idx = dirj * 2 + diri;
			auto cit = binary_term.insert(std::make_pair(std::make_pair(ccidj, ccidi), std::vector<double>(4)));
			if (cit.second == true) {
				cit.first->second[0] = 0.0;
				cit.first->second[1] = 0.0;
				cit.first->second[2] = 0.0;
				cit.first->second[3] = 0.0;
			}
			cit.first->second[idx] += sc;
		}

	}

	MRF_solver mrf(_gh.node_list.size());

	auto get_unary_score = [&](int i, int pos_neg) {
		return unary_term[2 * i + pos_neg];
	};

	for (int i = 0; i < _gh.node_list.size(); i++) {
		double val0 = get_unary_score(i, 0);
		double val1 = get_unary_score(i, 1);

		auto it = node_induced_dir.find(i);
		if (it != node_induced_dir.end()) {
			if (it->second) {
				val0 += 100000;
			}
			else {
				val1 += 100000;
			}
		}

		mrf.add_nodes(i, -1 * val0, -1 * val1);
	}
	for (auto it = binary_term.begin(); it != binary_term.end(); it++) {
		int cci = it->first.first;
		int ccj = it->first.second;
		double c00, c01, c10, c11;
		c00 = it->second[0];
		c01 = it->second[1];
		c10 = it->second[2];
		c11 = it->second[3];
		mrf.add_edge(cci, ccj, -1 * c00, -1 * c01, -1 * c10, -1 * c11);
	}

	ccdir.clear();
	mrf.solve(ccdir);

}

int Optimizer::get_max_feasiable_confid(std::vector<vector<double>>& sol, std::vector<int>& ccdir, int edgeid)
{

	std::vector<std::pair<int, double>> confid_poss;
	int degree = _gh.edge_list[_gh.non_manifold_curves[edgeid][0]].incident_components.size();
	int maxid = -1;
	double maxposs = -1.0;
	for (int i = 0; i < configuration_list_all_degrees[degree - 2].size(); i++) {
		if (sol[edgeid][i] < maxposs)continue;

		auto duallist = configuration_list_all_degrees[degree - 2][i];

		bool matched_flag = true;

		for (int j = 0; j < duallist.size(); j++) {
			auto dualedge = dual_edge_list_all_degrees[degree - 2][duallist[j]];
			int inci = dualedge.fi;
			int incj = dualedge.fj;
			int dir = dualedge.dir;

			int ccidi = _gh.edge_list[_gh.non_manifold_curves[edgeid][0]].incident_components[inci];


			if (dir != ccdir[ccidi]) {
				matched_flag = false;
				break;//ccidi's dir does not match ccdir	
			}

			if (incj == -1) {
				//dir == ccdir[ccidi]
				continue;
			};

			int ccidj = _gh.edge_list[_gh.non_manifold_curves[edgeid][0]].incident_components[incj];
			auto eh = _gh.edge_list[_gh.non_manifold_curves[edgeid][0]];
			int fidfi = _msh->GetEdges()[eh.original_edge_id].incident_faces[inci];
			int fidfj = _msh->GetEdges()[eh.original_edge_id].incident_faces[incj];

			BlackMesh::BlackMesh<double>::Triangle tr1 = _msh->GetTriangles()[fidfi];
			BlackMesh::BlackMesh<double>::Triangle tr2 = _msh->GetTriangles()[fidfj];
			int cohe_flg = if_face_coherant_based_on_edge(_msh->GetEdges()[eh.original_edge_id].vertices.first,
				_msh->GetEdges()[eh.original_edge_id].vertices.second,
				tr1, tr2);

			int induced_dirj = -1;

			if (dir == 0 && cohe_flg == 1) {
				induced_dirj = 0;
			}
			else if (dir == 1 && cohe_flg == 1) {
				induced_dirj = 1;
			}
			else if (dir == 0 && cohe_flg == 0) {
				induced_dirj = 1;
			}
			else if (dir == 1 && cohe_flg == 0) {
				induced_dirj = 0;
			}

#ifdef _OUTPUT_DEBUG_LOGS
			if (induced_dirj == -1) {
				std::cout << "Weird\n";
#if defined(_RECORD_ERROR_LOGS)
				(*_log_off) << _file_under_processing << " ERROR get_max_feasiable_confid\n";
				(*_log_off).flush();
#endif
			}
#endif
			if (induced_dirj != ccdir[ccidj]) {
				matched_flag = false;
				break;//ccidj'd dir is not matched
			}
		}

		if (matched_flag) {
			maxposs = sol[edgeid][i];
			maxid = i;
		}
	}
	if (maxid == -1) {
		std::cout << "Weird\n";
#if defined(_RECORD_ERROR_LOGS)
		(*_log_off) << _file_under_processing << " ERROR get_max_feasiable_confid\n";
		(*_log_off).flush();
#endif
	}
	return maxid;
}

void Optimizer::MRF_rounding(std::vector<vector<double>>& sol, std::map<int, int> *node, std::vector<vector<double>>& rounded)
{
	rounded.clear();
	std::vector<int> ccdir;
	solve_ccdir_on_current_sol(sol, node, ccdir);
	for (int i = 0; i < _gh.non_manifold_curves.size(); i++) {
		int id = get_max_feasiable_confid(sol, ccdir, i);

		rounded.push_back(vector<double>(sol[i].size(), 0));
		rounded[i][id] = 1;
	}
}


void Optimizer::update_curve_length()
{
	auto get_edge_len = [&](int nedgeid) {
		int edgeid = _gh.edge_list[nedgeid].original_edge_id;
		int vid1 = _msh->GetEdges()[edgeid].vertices.first;
		int vid2 = _msh->GetEdges()[edgeid].vertices.second;

		auto vp1 = _msh->GetVertices()[vid1].pos;
		auto vp2 = _msh->GetVertices()[vid2].pos;

		return sqrt(pow(vp1[0] - vp2[0], 2)
			+ pow(vp1[1] - vp2[1], 2)
			+ pow(vp1[2] - vp2[2], 2));
	};

	non_manifold_curve_length.clear();
	for (int i = 0; i < _gh.non_manifold_curves.size(); i++) {

		double curvelen = 0.0;
		for (int j = 0; j < _gh.non_manifold_curves[i].size(); j++) {
			auto edgeid = _gh.non_manifold_curves[i][j];
			curvelen += get_edge_len(edgeid);
		}

		non_manifold_curve_length.push_back(curvelen);
	}
}

void Optimizer::init()
{
	generate_transformation_matrices();
}

//stack edges to curves
void Optimizer::stack_non_manifold_curves(std::vector<std::vector<int>> &nonmanifold_curves) {
	auto build_oriedge_to_graphedge_map = [](map<int, int> &res, BlackMesh::BlackMesh<double> *msh, Primal_Dual_graph *gh) {
		res.clear();

		for (int i = 0; i < gh->edge_list.size(); i++) {

			res.insert(std::make_pair(gh->edge_list[i].original_edge_id, i));

		}

#ifdef _OUTPUT_DEBUG_LOGS
		for (int i = 0; i < gh->edge_list.size(); i++) {
			if (res[gh->edge_list[i].original_edge_id] != i) {
				std::cout << "The map between edges is not matched\n";
#if defined(_RECORD_ERROR_LOGS)
				(*_log_off) << _file_under_processing << " The map between edges is not matched\n";
				(*_log_off).flush();
#endif
			}
		}
#endif
	};

	auto get_incident_nedge_list = [](int edgeid_msh, std::vector<int> &list, std::map<int, int> &old_to_new, BlackMesh::BlackMesh<double> *msh) {
		//get all incident nedges connecting vtxid1 vtxid2
		list.clear();
		//two vertices
		int vtxid1 = msh->GetEdges()[edgeid_msh].vertices.first;
		int vtxid2 = msh->GetEdges()[edgeid_msh].vertices.second;

		int counter;

		counter = 0;
		for (int i = 0; i < msh->GetVertices()[vtxid1].incident_edges.size(); i++) {
			int edgeid = msh->GetVertices()[vtxid1].incident_edges[i];
			if (edgeid == edgeid_msh)continue;
			std::map<int, int>::iterator it = old_to_new.find(edgeid);

			if (it != old_to_new.end()) {

				list.push_back(it->second);
				counter++;
			}
#ifdef _OUTPUT_DEBUG_LOGS
			else {
				if (msh->GetEdges()[edgeid].is_nonmanifold_edge == 1) {
					std::cout << "Query edge do not exist in the mapping\n";
#if defined(_RECORD_ERROR_LOGS)
					(*_log_off) << _file_under_processing << " ERROR Query edge do not exist in the mapping\n";
					(*_log_off).flush();
#endif
				}
			}
#endif
		}
		//if (counter > 1)return false;

		counter = 0;
		for (int i = 0; i < msh->GetVertices()[vtxid2].incident_edges.size(); i++) {
			int edgeid = msh->GetVertices()[vtxid2].incident_edges[i];
			if (edgeid == edgeid_msh)continue;
			std::map<int, int>::iterator it = old_to_new.find(edgeid);

			if (it != old_to_new.end()) {

				list.push_back(it->second);
				counter++;
			}
#ifdef _OUTPUT_DEBUG_LOGS
			else {
				if (msh->GetEdges()[edgeid].is_nonmanifold_edge == 1) {
					std::cout << "Query edge do not exist in the mapping\n";
#if defined(_RECORD_ERROR_LOGS)
					(*_log_off) << _file_under_processing << " ERROR Query edge do not exist in the mapping\n";
					(*_log_off).flush();
#endif
				}
			}
#endif
		}
		//if (counter > 1)return false;

		//return true;
	};

	auto if_edge_have_same_incident_cc_number = [](int id1, int id2, BlackMesh::BlackMesh<double> *msh) {
		return msh->GetEdges()[id1].incident_faces.size() == msh->GetEdges()[id2].incident_faces.size();
	};

	auto if_two_face_identical_or_coherant = [&](int id1, int id2, BlackMesh::BlackMesh<double> *msh) {
		if (id1 == id2) return true;

		std::vector<int> same_vtx_id;
		for (int i = 0; i < msh->GetTriangles()[id1].vertices.size(); i++) {
			for (int j = 0; j < msh->GetTriangles()[id2].vertices.size(); j++) {
				if (msh->GetTriangles()[id1].vertices[i] == j < msh->GetTriangles()[id2].vertices[j]) {
					same_vtx_id.push_back(msh->GetTriangles()[id1].vertices[i]);
					break;
				}
			}
		}

		if (same_vtx_id.size() < 2) return false;
		if (same_vtx_id.size() > 2) {
#ifdef _OUTPUT_DEBUG_LOGS
			std::cout << "Weird things happen\n";
#if defined(_RECORD_ERROR_LOGS)
			(*_log_off) << _file_under_processing << " ERROR Weird things happen\n";
			(*_log_off).flush();
#endif
#endif
		}

		return 1 == if_face_coherant_based_on_edge(same_vtx_id[0], same_vtx_id[1], msh->GetTriangles()[id1], msh->GetTriangles()[id2]);
	};

	auto if_edge_incident_cc_are_same_and_coherant = [&](int ghidQuery, int ghidBase, BlackMesh::BlackMesh<double> *msh) {
		int id1 = _gh.edge_list[ghidQuery].original_edge_id;
		int id2 = _gh.edge_list[ghidBase].original_edge_id;



		//
		if (if_edge_have_same_incident_cc_number(id1, id2, msh) == false)return false;

		//if the set of incident id are same
		for (int i = 0; i < msh->GetEdges()[id1].incident_faces.size(); i++) {
			int ccidi = msh->GetTriangles()[msh->GetEdges()[id1].incident_faces[i]].component_id;
			bool matched = false;
			for (int j = 0; j < msh->GetEdges()[id2].incident_faces.size(); j++) {
				int ccidj = msh->GetTriangles()[msh->GetEdges()[id2].incident_faces[j]].component_id;
				if (ccidi == ccidj) {
					//should take the case where two face have same ccid into consideration
					matched = true;
					break;

				}
			}
			if (matched == false) return false;
		}


		//if coherent, if so, change the storage in mesh and gh
		//get share vtx
		int share_vtx = -1;
		if (_msh->GetEdges()[id1].vertices.first == _msh->GetEdges()[id2].vertices.first) {
			share_vtx = _msh->GetEdges()[id1].vertices.first;
		}
		else if (_msh->GetEdges()[id1].vertices.first == _msh->GetEdges()[id2].vertices.second) {
			share_vtx = _msh->GetEdges()[id1].vertices.first;
		}
		else if (_msh->GetEdges()[id1].vertices.second == _msh->GetEdges()[id2].vertices.first) {
			share_vtx = _msh->GetEdges()[id1].vertices.second;
		}
		else if (_msh->GetEdges()[id1].vertices.second == _msh->GetEdges()[id2].vertices.second) {
			share_vtx = _msh->GetEdges()[id1].vertices.second;
		}

		//
		std::vector<std::vector<int>> eq_face_set;
		eq_face_set.resize(_gh.edge_list[ghidBase].incident_components.size());
		for (int i = 0; i < eq_face_set.size(); i++) {
			int start_face = msh->GetEdges()[id2].incident_faces[i];
			std::queue<int> list;
			list.push(start_face);
			eq_face_set[i].push_back(start_face);
			std::set<int> visited_face;
			visited_face.insert(start_face);

			while (list.size() != 0) {
				int this_face_id = list.front();
				list.pop();

				//std::vector<int> incident_faces_to_this_face;
				auto fh = msh->GetTriangles()[this_face_id];

				for (int j = 0; j < fh.incident_edges.size(); j++) {
					int edgeid = fh.incident_edges[j];
					if (edgeid == id1 || edgeid == id2)continue;

					for (int m = 0; m < msh->GetEdges()[edgeid].incident_faces.size(); m++) {
						int neiFace = msh->GetEdges()[edgeid].incident_faces[m];

						//if neiFace have share_vtx as vtx
						bool have_share_vtx = false;
						for (int mm = 0; mm < msh->GetTriangles()[neiFace].vertices.size(); mm++) {
							if (msh->GetTriangles()[neiFace].vertices[mm] == share_vtx) {
								have_share_vtx = true;
								break;
							}
						}
						if (have_share_vtx == false) continue;

						auto it = visited_face.insert(neiFace);
						if (it.second == true) {
							//do not exist
							eq_face_set[i].push_back(neiFace);
							list.push(neiFace);
						}
					}
				}
			}

		}

		//if all incident faces of queryid could be found in eqlist, return true;
		std::map<int, int> newIncidentId;//for each query id, which incident id in base edge is equal
		bool all_matched = true;
		for (int i = 0; i < msh->GetEdges()[id1].incident_faces.size(); i++) {
			int query_face_id = msh->GetEdges()[id1].incident_faces[i];
			bool found_flag = false;
			for (int classid = 0; classid < eq_face_set.size(); classid++) {
				for (int eleid = 0; eleid < eq_face_set[classid].size(); eleid++) {
					if (query_face_id == eq_face_set[classid][eleid]) {
						found_flag = true;
						newIncidentId.insert(std::make_pair(classid, i));
						break;
					}
				}
				if (found_flag)break;
			}

			if (found_flag == false) {
				all_matched = false;
				break;
			}
		}

		if (all_matched) {
			auto ori_gh_list = _gh.edge_list[ghidQuery].incident_components;
			auto ori_msh_list = msh->GetEdges()[id1].incident_faces;

			for (int i = 0; i < ori_gh_list.size(); i++) {
				_gh.edge_list[ghidQuery].incident_components[i] = ori_gh_list[newIncidentId[i]];
				(msh->GetEdgesEditable()[id1].incident_faces)[i] = ori_msh_list[newIncidentId[i]];
			}
			return true;
		}
		else {
			return false;
		}


	};


	auto if_two_edge_share_vtx_have_more_than_two_nedge = [&](int id1, int id2) {
		int share_vtx = -1;
		if (_msh->GetEdges()[id1].vertices.first == _msh->GetEdges()[id2].vertices.first) {
			share_vtx = _msh->GetEdges()[id1].vertices.first;
		}
		else if (_msh->GetEdges()[id1].vertices.first == _msh->GetEdges()[id2].vertices.second) {
			share_vtx = _msh->GetEdges()[id1].vertices.first;
		}
		else if (_msh->GetEdges()[id1].vertices.second == _msh->GetEdges()[id2].vertices.first) {
			share_vtx = _msh->GetEdges()[id1].vertices.second;
		}
		else if (_msh->GetEdges()[id1].vertices.second == _msh->GetEdges()[id2].vertices.second) {
			share_vtx = _msh->GetEdges()[id1].vertices.second;
		}

		int ncounter = 0;
		for (int s = 0; s < _msh->GetVertices()[share_vtx].incident_edges.size(); s++) {
			auto flg = _msh->GetEdges()[_msh->GetVertices()[share_vtx].incident_edges[s]].is_nonmanifold_edge;
			if (flg == true)ncounter++;
		}

		return ncounter != 2;
	};

	auto flood_pass_approved = [&](int base_id_gh, std::vector<int> &input, std::vector<int> &output,
		BlackMesh::BlackMesh<double> *msh, Primal_Dual_graph *gh) {

		output.clear();
		int base_id_msh = gh->edge_list[base_id_gh].original_edge_id;
		for (int i = 0; i < input.size(); i++) {
			int query_id_gh = input[i];



			if (query_id_gh == base_id_gh)continue;

			int query_id_msh = gh->edge_list[query_id_gh].original_edge_id;

			bool ignore_flag = if_two_edge_share_vtx_have_more_than_two_nedge(query_id_msh, base_id_msh);
			if (ignore_flag)continue;

			bool flg = if_edge_incident_cc_are_same_and_coherant(query_id_gh, base_id_gh, msh);

			if (flg) output.push_back(query_id_gh);
		}
	};

	//std::vector<std::vector<int>> nonmanifold_curves;
	std::vector<bool> if_edge_moved;
	if_edge_moved.resize(_gh.edge_list.size(), false);

	std::map<int, int> old_to_new;
	build_oriedge_to_graphedge_map(old_to_new, _msh, &_gh);

	for (int i = 0; i < if_edge_moved.size(); i++) {

		if (if_edge_moved[i]) continue;

		std::vector<int> this_curve;

		std::queue<int> list;
		list.push(i);
		if_edge_moved[i] = true;
		this_curve.push_back(i);

		while (list.size() != 0) {
			int current_id_gh = list.front();
			list.pop();

			std::vector<int> incident_nedge_list;
			std::vector<int> approved_nedge_list;
			get_incident_nedge_list(_gh.edge_list[current_id_gh].original_edge_id, incident_nedge_list, old_to_new, _msh);



			flood_pass_approved(current_id_gh, incident_nedge_list, approved_nedge_list, _msh, &_gh);

			for (int i = 0; i < approved_nedge_list.size(); i++) {
				int idd = approved_nedge_list[i];
				if (if_edge_moved[idd] == false) {
					list.push(idd);
					if_edge_moved[idd] = true;
					this_curve.push_back(idd);
				}
			}



		}
		nonmanifold_curves.push_back(this_curve);
	}

#ifdef _OUTPUT_MERGED_NONMANIFOLD_EDGES
	for (int i = 0; i < nonmanifold_curves.size(); i++) {
		char buffer[128];
		sprintf(buffer, "Curv\\_debug_nedge_curve_pt%d.obj", i);
		ofstream off(buffer);
		for (int s = 0; s < _msh->GetNumVertices(); s++) {
			off << "v " << _msh->GetVertices()[s].pos[0] << " " <<
				_msh->GetVertices()[s].pos[1] << " " <<
				_msh->GetVertices()[s].pos[2] << " \n";
		}

		for (int l = 0; l < nonmanifold_curves[i].size(); l++) {
			int ide_msh = _gh.edge_list[nonmanifold_curves[i][l]].original_edge_id;
			off << "l " << _msh->GetEdges()[ide_msh].vertices.first + 1 << " " << _msh->GetEdges()[ide_msh].vertices.second + 1 << "\n";
		}
		off.close();
	}
#endif
}


//assemble
void Optimizer::assemble_orientation_cons(int edgeid)
{
	Primal_Dual_graph::Nonmanifold_edge eh = _gh.edge_list[edgeid];

	int degree = eh.incident_components.size();

	std::vector<int> idxrow_conf2face, idxcol_conf2face;
	for (int m = 0; m < configuration_list_all_degrees[degree - 2].size(); m++) {
		std::vector<int> conf = configuration_list_all_degrees[degree - 2][m];


		for (int s = 0; s < conf.size(); s++) {
			//each dual edge in this conf
			int idfi = dual_edge_list_all_degrees[degree - 2][conf[s]].fi;
			int idfj = dual_edge_list_all_degrees[degree - 2][conf[s]].fj;
			int dir = dual_edge_list_all_degrees[degree - 2][conf[s]].dir;//0-positive
																		  //assert(idfi < idfj);

			if (idfi != -1 && idfj != -1) {
				//original face id
				int fidfi = _msh->GetEdges()[eh.original_edge_id].incident_faces[idfi];
				int fidfj = _msh->GetEdges()[eh.original_edge_id].incident_faces[idfj];

				BlackMesh::BlackMesh<double>::Triangle tr1 = _msh->GetTriangles()[fidfi];
				BlackMesh::BlackMesh<double>::Triangle tr2 = _msh->GetTriangles()[fidfj];

				int flg = if_face_coherant_based_on_edge(_msh->GetEdges()[eh.original_edge_id].vertices.first,
					_msh->GetEdges()[eh.original_edge_id].vertices.second,
					tr1, tr2);

				assert(flg != -1);

				std::map<int, std::vector<orientation_Expr_with_id>> tmp;
				ori_cons_with_confid[eh.incident_components[idfi]].insert(std::make_pair(edgeid, tmp));
				ori_cons_with_confid[eh.incident_components[idfj]].insert(std::make_pair(edgeid, tmp));

				ori_cons_with_confid[eh.incident_components[idfi]][edgeid].insert(std::make_pair(idfi, std::vector<orientation_Expr_with_id>{}));
				ori_cons_with_confid[eh.incident_components[idfj]][edgeid].insert(std::make_pair(idfj, std::vector<orientation_Expr_with_id>{}));

				if (dir == 0) {
					ori_cons_with_confid[eh.incident_components[idfi]][edgeid][idfi].push_back(orientation_Expr_with_id(true, edgeid, m));
					if (flg == 1) {
						ori_cons_with_confid[eh.incident_components[idfj]][edgeid][idfj].push_back(orientation_Expr_with_id(true, edgeid, m));
					}
					else {
						ori_cons_with_confid[eh.incident_components[idfj]][edgeid][idfj].push_back(orientation_Expr_with_id(false, edgeid, m));
					}
				}
				else {
					ori_cons_with_confid[eh.incident_components[idfi]][edgeid][idfi].push_back(orientation_Expr_with_id(false, edgeid, m));
					if (flg == 1) {
						ori_cons_with_confid[eh.incident_components[idfj]][edgeid][idfj].push_back(orientation_Expr_with_id(false, edgeid, m));
					}
					else {
						ori_cons_with_confid[eh.incident_components[idfj]][edgeid][idfj].push_back(orientation_Expr_with_id(true, edgeid, m));
					}
				}


			}
			else {
				int idf = idfi == -1 ? idfj : idfi;
				if (dir == 0) {
					ori_cons_with_confid[eh.incident_components[idf]][edgeid][idf].push_back(orientation_Expr_with_id(true, edgeid, m));
				}
				else {
					ori_cons_with_confid[eh.incident_components[idf]][edgeid][idf].push_back(orientation_Expr_with_id(false, edgeid, m));
				}
			}
		}
	}
}

void Optimizer::assemble_orientation_cons_curve_based(int curveid)
{
	int edgeid = _gh.non_manifold_curves[curveid][0];//use the first edge of each curve as present
	Primal_Dual_graph::Nonmanifold_edge eh = _gh.edge_list[edgeid];

	int degree = eh.incident_components.size();

	std::vector<int> idxrow_conf2face, idxcol_conf2face;
	for (int m = 0; m < configuration_list_all_degrees[degree - 2].size(); m++) {
		std::vector<int> conf = configuration_list_all_degrees[degree - 2][m];


		for (int s = 0; s < conf.size(); s++) {
			//each dual edge in this conf
			int dualedgeid = conf[s];

			int idfi = dual_edge_list_all_degrees[degree - 2][conf[s]].fi;
			int idfj = dual_edge_list_all_degrees[degree - 2][conf[s]].fj;
			int dir = dual_edge_list_all_degrees[degree - 2][conf[s]].dir;//0-positive
																		  //assert(idfi < idfj);

			if (idfi != -1 && idfj != -1) {
				//original face id
				int fidfi = _msh->GetEdges()[eh.original_edge_id].incident_faces[idfi];
				int fidfj = _msh->GetEdges()[eh.original_edge_id].incident_faces[idfj];

				int ccidi = eh.incident_components[idfi];
				int ccidj = eh.incident_components[idfj];

				//int minccid = min(ccidi, ccidj);
				//int maxccid = max(ccidi, ccidj);
#ifdef _OUTPUT_DEBUG_LOGS
				if (ccidi != _msh->GetTriangles()[fidfi].component_id) {
					std::cout << "Weird\n";
#if defined(_RECORD_ERROR_LOGS)
					(*_log_off) << _file_under_processing << " ERROR Weird in orientation assemble\n";
					(*_log_off).flush();
#endif
				}
				if (ccidj != _msh->GetTriangles()[fidfj].component_id) {
					std::cout << "Weird\n";
#if defined(_RECORD_ERROR_LOGS)
					(*_log_off) << _file_under_processing << " ERROR Weird in orientation assemble\n";
					(*_log_off).flush();
#endif
				}
#endif 

				BlackMesh::BlackMesh<double>::Triangle tr1 = _msh->GetTriangles()[fidfi];
				BlackMesh::BlackMesh<double>::Triangle tr2 = _msh->GetTriangles()[fidfj];

				int flg = if_face_coherant_based_on_edge(_msh->GetEdges()[eh.original_edge_id].vertices.first,
					_msh->GetEdges()[eh.original_edge_id].vertices.second,
					tr1, tr2);

				assert(flg != -1);

				std::map<int, std::vector<orientation_Expr_with_id>> tmp;
				ori_cons_with_confid[eh.incident_components[idfi]].insert(std::make_pair(curveid, tmp));
				ori_cons_with_confid[eh.incident_components[idfj]].insert(std::make_pair(curveid, tmp));

				ori_cons_with_confid[eh.incident_components[idfi]][curveid].insert(std::make_pair(idfi, std::vector<orientation_Expr_with_id>{}));
				ori_cons_with_confid[eh.incident_components[idfj]][curveid].insert(std::make_pair(idfj, std::vector<orientation_Expr_with_id>{}));

				if (dir == 0) {
					ori_cons_with_confid[eh.incident_components[idfi]][curveid][idfi].push_back(orientation_Expr_with_id(true, curveid, m));
					if (flg == 1) {
						ori_cons_with_confid[eh.incident_components[idfj]][curveid][idfj].push_back(orientation_Expr_with_id(true, curveid, m));
						if (ccidi != ccidj) {
							insert_to_binary_score(curveid, m, dualedgeid, idfi, idfj, 0, 0);
							insert_to_unary_score(curveid, m, ccidi, 0);
							insert_to_unary_score(curveid, m, ccidj, 0);

						}
						else {
							insert_to_unary_from_binary_score(curveid, m, dualedgeid, idfi, idfj, 0);
							insert_to_unary_score(curveid, m, ccidi, 0);

						}
					}
					else {
						ori_cons_with_confid[eh.incident_components[idfj]][curveid][idfj].push_back(orientation_Expr_with_id(false, curveid, m));
						if (ccidi != ccidj) {
							insert_to_binary_score(curveid, m, dualedgeid, idfi, idfj, 0, 1);
							insert_to_unary_score(curveid, m, ccidi, 0);
							insert_to_unary_score(curveid, m, ccidj, 1);
						}
						else {

						}
					}
				}
				else {
					ori_cons_with_confid[eh.incident_components[idfi]][curveid][idfi].push_back(orientation_Expr_with_id(false, curveid, m));
					if (flg == 1) {
						ori_cons_with_confid[eh.incident_components[idfj]][curveid][idfj].push_back(orientation_Expr_with_id(false, curveid, m));
						if (ccidi != ccidj) {
							insert_to_binary_score(curveid, m, dualedgeid, idfi, idfj, 1, 1);
							insert_to_unary_score(curveid, m, ccidi, 1);
							insert_to_unary_score(curveid, m, ccidj, 1);
						}
						else {
							insert_to_unary_from_binary_score(curveid, m, dualedgeid, idfi, idfj, 1);
							insert_to_unary_score(curveid, m, ccidi, 1);
						}
					}
					else {
						ori_cons_with_confid[eh.incident_components[idfj]][curveid][idfj].push_back(orientation_Expr_with_id(true, curveid, m));
						if (ccidi != ccidj) {
							insert_to_binary_score(curveid, m, dualedgeid, idfi, idfj, 1, 0);
							insert_to_unary_score(curveid, m, ccidi, 1);
							insert_to_unary_score(curveid, m, ccidj, 0);
						}
						else {
						}
					}
				}


			}
			else {
				int idf = idfi == -1 ? idfj : idfi;
				if (dir == 0) {
					ori_cons_with_confid[eh.incident_components[idf]][curveid][idf].push_back(orientation_Expr_with_id(true, curveid, m));
					insert_to_unary_score(curveid, m, eh.incident_components[idf], 0);
				}
				else {
					ori_cons_with_confid[eh.incident_components[idf]][curveid][idf].push_back(orientation_Expr_with_id(false, curveid, m));
					insert_to_unary_score(curveid, m, eh.incident_components[idf], 1);
				}
			}
		}
	}
}

void Optimizer::assemble_contraints_on_one_nonmanifold_edge(int id)
{
	//append orientation cons
	assemble_orientation_cons_curve_based(id);

}

void Optimizer::assemble_contraints() {
	//assemble all global connectivity cons and orientation cons

	ori_cons_with_confid.clear();
	ori_cons_with_confid.resize(_gh.node_list.size());

	//also assemble
	binary_score.clear();
	unary_score_from_binary.clear();
	unary_score.clear();
	unary_score.resize(_gh.node_list.size(), Unary_Element(_gh.node_list.size()));


	for (int i = 0; i < _gh.non_manifold_curves.size(); i++) {
		assemble_contraints_on_one_nonmanifold_edge(i);
	}

}

void Optimizer::assemble_imcompatible_confs()
{
	imcompatible_confs.clear();
	
	std::vector<bool> edge_check_list;
	edge_check_list.resize(_gh.non_manifold_curves.size(), false);

	for (int i = 0; i < _gh.non_manifold_curves.size(); i++) {
		int edgeid = _gh.non_manifold_curves[i][0];

		auto incident_list = _gh.edge_list[edgeid].incident_components;

		std::set<int> incident_group(incident_list.begin(), incident_list.end());

		int deg = incident_list.size();

		if (incident_group.size() != incident_list.size()) {
			//need check
			edge_check_list[i] = true;

		}

	}

	//	std::vector<map<int, map<int,std::vector<orientation_Expr_with_id>>>> 
	//ori_cons_with_confid;
	//ccid,edge id, incident id, cons
	
	//edgeid,conf-ccid +/-
	std::map<std::pair<int, int>, std::map<int, bool>> conf_cc_dir;

	for (int ccid=0; ccid!= ori_cons_with_confid.size(); ccid++) {
		for (auto it = ori_cons_with_confid[ccid].begin(); 
			it != ori_cons_with_confid[ccid].end(); it++) {
		
			int curvid = it->first;
			if (edge_check_list[curvid]) {
				for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
					//each incident id
					for (auto ittt = itt->second.begin(); ittt != itt->second.end(); ittt++) {
						int confid = ittt->confid;
						auto fit=conf_cc_dir.insert(std::make_pair(std::make_pair(curvid, confid), std::map<int, bool>{}));

						if (ittt->positive_flg) {
							auto fitt = fit.first->second.insert(std::make_pair(ccid, true));
							if (!fitt.second) {
								if (fitt.first->second != true) {
									imcompatible_confs.push_back(std::make_pair(curvid, confid));
								}
							}
						}
						else {
							auto fitt = fit.first->second.insert(std::make_pair(ccid, false));
							if (!fitt.second) {
								if (fitt.first->second != false) {
									imcompatible_confs.push_back(std::make_pair(curvid, confid));
								}
							}
						}
					}
				}
			}
		}
	}
}

void Optimizer::append_imcp_confs(map<std::pair<int, int>, int>& mp, Mosek_c_wrapper & solver, std::vector<std::pair<int, int>>& imcompatible_confs)
{


		std::vector<std::vector<int>> vid;
		std::vector<std::vector<double>> val;

		for (auto it = imcompatible_confs.begin(); it != imcompatible_confs.end(); it++) {
			int edgeid = it->first;
			int confid = it->second;
			vid.push_back(std::vector<int>{mp[std::make_pair(edgeid, confid)]});
			val.push_back(std::vector<double>{1.0});
		}

		if (vid.size() > 0)
			solver.append_fixed_var_cons(vid, val,0.0);
	
}

//append
void Optimizer::append_cons_on_all_edges(std::map<int, int>& node)
{

}


//solve
#define _ENTROPY_SIMIALR_EPS_ 1e-5
void Optimizer::branch_and_bound(std::vector<std::vector<double>> &xopt) {
	//the Energy will be MINIMIZED
	xopt.clear();

	//
	auto entropy_of_vector = [](const std::vector<double> &x) {

		double entropy = 0.0;
		for (int i = 0; i < x.size(); i++) {
			if (x[i] > 1e-5) {
				entropy += -1.0 * x[i]*log(x[i]);
			}
		}

		return entropy;
	};

	//initialize
	double upper_bound = 100000000;
	double lower_bound = -100000000;
	double opt_eneygy = -0;
	std::map<int, int> node;//edgeid, configuration

	std::queue<std::map<int, int>> list;
	list.push(node);


#ifdef _USE_TIMER
	clock_t t1, t2;
	t1 = clock();
	double END_TIMER_MIN = TIME_MIN * 60.0 * (double)(CLOCKS_PER_SEC);
	double END_TIMER_MAX = TIME_MAX * 60.0 * (double)(CLOCKS_PER_SEC);
#endif

	int iter_counter = 0;
	clock_t t_opt;
	int iter_counter_opt = 0;


	int opt_reached_counter = 0;


	while (list.size() != 0) {
		//(node, G0)=pop(list)
		//int ls_size = list.size();
		node = list.front();
		list.pop();

		iter_counter++;


		//solve solution on the graph induced by node
		std::vector<std::vector<double>> solutions;
		auto eng = solve_problem_with_fixed_solution_capi(node, solutions);

		if (eng[0] == -1 && eng[4] == -1) {
			//solver failed
#if defined(_RECORD_ERROR_LOGS)
			(*_log_off) << _file_under_processing << " ERROR solver failed in solving\n";
			(*_log_off).flush();
#endif
			return;
		}

		lower_bound = -1 * eng[0];
		//Lower_bound=E(x)
#ifdef _OUTPUT_BRANCH_BOUND_LOG
		std::cout << "**Solved, Lowerbound is " << lower_bound <<
			", ConnBd is " << -eng[1] << ", VisSc is " << -eng[2] <<
			", Uterm " << -eng[3] << ", Bterm " << -eng[4] << "\n";
#endif

		std::vector<std::vector<double>> solutions_rounded;

		MRF_rounding(solutions, &node, solutions_rounded);

		//
		if (if_x_feasible(solutions_rounded)) {
#ifdef _OUTPUT_BRANCH_BOUND_LOG
			std::cout << "****Rounded solution is feasible \n";
#endif
			auto engr = algebraic_conn_bound_capi(solutions_rounded);

			if (engr[0] == -1 && engr[4] == -1) {
				//solver failed
#if defined(_RECORD_ERROR_LOGS)
				(*_log_off) << _file_under_processing << " ERROR solver failed in rounding energy evaluation\n";
				(*_log_off).flush();
#endif

				return;
			}

			double rounded_energy = -1 * engr[0];
#ifdef _OUTPUT_BRANCH_BOUND_LOG
			std::cout << "****Rounded Energy is: " << rounded_energy <<
				", ConnBd is " << -engr[1] << ", VisSc is " << -engr[2] <<
				", Uterm " << -engr[3] << ", Bterm " << -engr[4] << "\n";

			if (lower_bound - rounded_energy > 1e-5) {
#if defined(_RECORD_ERROR_LOGS) 
				(*_log_off) << _file_under_processing <<
					" ERROR rounded energy become better o/r " << eng[0] << "/" << engr[0] << "\n";
				(*_log_off).flush();
#endif
			}
#endif
			if (rounded_energy < upper_bound) {
				//upper_bound = -1 * eng[0];//
				upper_bound = rounded_energy;
				xopt = solutions_rounded;
				opt_eneygy = rounded_energy;
				t_opt = clock();
				iter_counter_opt = iter_counter;
#ifdef _ITERATE_ONLY_ONE
				break;
#endif

				if (_opt_max_iter != -1 && _opt_max_iter <= ++opt_reached_counter) {
					break;
				}

#ifdef _OUTPUT_BRANCH_BOUND_LOG
				std::cout << "******Upperbound updated, " << upper_bound << "\n";

#endif
			}
		}
		else {

			std::cout << "****Rounded solution is NOT feasible, SHOULD NOT HAPPEN HERE!!! \n";
#if defined(_RECORD_ERROR_LOGS) 
			(*_log_off) << _file_under_processing << " ERROR Rounded solution is NOT feasible\n";
			(*_log_off).flush();
#endif

		}

		if (lower_bound - upper_bound > 1e-7) {
			//do nothing
#ifdef _OUTPUT_BRANCH_BOUND_LOG
			std::cout << "****Lower bound larger than UpperBound, do nothing\n";
#endif
		}
		else {
			//find x with maximal entropy
			double max_entropy = -100000000;
			int max_id = -1;
			for (int s = 0; s < solutions.size(); s++) {
				if (node.find(s) != node.end())continue;
				double ent = entropy_of_vector(solutions[s]);
				if (ent > max_entropy) {
					max_id = s;
					max_entropy = ent;
				}
			}
#ifdef _OUTPUT_BRANCH_BOUND_LOG
			std::cout << "****Element " << max_id << " will be branched\n";
#endif
			if (max_id != -1) {
				std::pair<std::map<int, int>::iterator, bool> ret = node.insert(std::make_pair(max_id, -1));
				if (ret.second == false || solutions[max_id].size() == 0) {
					std::cout << "******Already Branched\n";
				}
				else {
					for (int j = 0; j < solutions[max_id].size(); j++) {
						node[max_id] = j;
						bool pushflg = if_branched_solution_meet_cons(node);
						//bool pushflg = true;
						if (pushflg) {
							list.push(node);
						}
						else {
							//	std::cout << "******To-branched "<<j<<" already violate OC\n";
						}
					}
				}
			}

		}

#ifdef _USE_TIMER
		t2 = clock();
		if (t2 - t1 > END_TIMER_MIN&&xopt.size() > 0) {
			std::cout << "Max Time Limit Reached\n";
			std::cout << "Energy " << opt_eneygy << "\n";
			std::cout << "Opt reached at iteration " << iter_counter_opt
				<< " ,time " << (t_opt - t1) / (60.0 * (double)(CLOCKS_PER_SEC)) << "\n";
			return;
		}
		else if (t2 - t1 > END_TIMER_MAX) {
			std::cout << "Energy " << opt_eneygy << "\n";
			std::cout << "Max Time Limit Reached\n";
			std::cout << "Opt reached at iteration " << iter_counter_opt
				<< " ,time " << (t_opt - t1) / (60.0 * (double)(CLOCKS_PER_SEC)) << "\n";
			return;
		}

#endif
	}

}

double  Optimizer::solve_problem_with_fixed_solution(std::map<int, int>& node, std::vector<std::vector<double>>& result_solution) {

	return 0.0;

}

////generated mosek wrapper needed data
//ccid -the lap index of each constraint
//vid-cons id, variable id
//val-cons id, coeff
void Optimizer::append_orientation_cons(map<std::pair<int, int>, int>& mp, Mosek_c_wrapper &solver)
{
	//ccid,edgeid,confid-x
	std::vector<std::vector<std::vector<int>>> ori_plus, ori_minus;
	for (int i = 0; i < ori_cons_with_confid.size(); i++) {
		std::vector<std::vector<int>> ori_plus_ci, ori_minus_ci;
		std::map<int, map<int, std::vector<orientation_Expr_with_id>>>::iterator it = ori_cons_with_confid[i].begin();
		for (; it != ori_cons_with_confid[i].end(); it++) {//go through all edges
			int edgeid = it->first;

			map<int, std::vector<orientation_Expr_with_id>>::iterator diff_incident_faces_it = it->second.begin();
			for (; diff_incident_faces_it != it->second.end(); diff_incident_faces_it++) {//each edge may have differnt incident faces as ci

				std::vector<int> ori_plus_ci_edgeit, ori_minus_ci_edgeit;
				bool pushed1 = false, pushed2 = false;

				//std::map<int, int>::iterator nit = node.find(edgeid);

					//not fixed
				for (int j = 0; j < diff_incident_faces_it->second.size(); j++) {
					int confid = (diff_incident_faces_it->second)[j].confid;
					bool dir = (diff_incident_faces_it->second)[j].positive_flg;


					if (dir) {
						//positive
						ori_plus_ci_edgeit.push_back(mp[std::make_pair(edgeid, confid)]);
						pushed1 = true;

					}
					else {
						//negative
						ori_minus_ci_edgeit.push_back(mp[std::make_pair(edgeid, confid)]);
						pushed2 = true;
					}
				}



				if (pushed1)
					ori_plus_ci.push_back(ori_plus_ci_edgeit);
				if (pushed2)
					ori_minus_ci.push_back(ori_minus_ci_edgeit);
			}


		}

		ori_plus.push_back(ori_plus_ci);
		ori_minus.push_back(ori_minus_ci);
	}

	//orientation constraints
	//if dual edge induced direction follow component normal, then the contraint should be pushed to positive slot; if not, negative
	//the dual edge positive rely on the direction of incident face with smaller incident idx in the incident list
	int cons_counter = 0;
	std::vector<map<int, double>> trilist;//use trilist to combine duplicate term.//rowid-colid-val
	for (int m = 0; m < ori_plus.size(); m++) {
		int siz = ori_plus[m].size();
		for (int i = 0; i < siz; i++) {
			int j = (i + 1) % siz;
			trilist.push_back(map<int, double>{});
			int consid = trilist.size();
			consid--;


			for (int ss = 0; ss < ori_plus[m][i].size(); ss++) {
				auto ret = trilist[consid].insert(std::make_pair(ori_plus[m][i][ss], 1.0));
				if (ret.second == false) {
					//already exist
					ret.first->second += 1.0;
				}

			}
			for (int ss = 0; ss < ori_plus[m][j].size(); ss++) {
				auto ret = trilist[consid].insert(std::make_pair(ori_plus[m][j][ss], -1.0));
				if (ret.second == false) {
					//already exist
					ret.first->second += -1.0;
				}
			}

		}
	}

	for (int m = 0; m < ori_minus.size(); m++) {
		int siz = ori_minus[m].size();
		for (int i = 0; i < siz; i++) {
			int j = (i + 1) % siz;
			trilist.push_back(map<int, double>{});
			int consid = trilist.size();
			consid--;



			for (int ss = 0; ss < ori_minus[m][i].size(); ss++) {
				auto ret = trilist[consid].insert(std::make_pair(ori_minus[m][i][ss], 1.0));
				if (ret.second == false) {
					//already exist
					ret.first->second += 1.0;
				}

			}
			for (int ss = 0; ss < ori_minus[m][j].size(); ss++) {
				auto ret = trilist[consid].insert(std::make_pair(ori_minus[m][j][ss], -1.0));
				if (ret.second == false) {
					//already exist
					ret.first->second += -1.0;
				}
			}

		}
	}

	std::vector<std::vector<int>> vid;
	std::vector<std::vector<double>> val;

	for (int i = 0; i < trilist.size(); i++) {
		vid.push_back(std::vector<int>{});
		val.push_back(std::vector<double>{});
		for (auto it = trilist[i].begin(); it != trilist[i].end(); it++) {
			if (abs(it->second - 0) < 1e-7) {
#ifndef _DISABLE_FILLING_WARNING
				cout << "problem,";
#endif
#if defined(_RECORD_ERROR_LOGS)
				//(*_log_off) << _file_under_processing << " ERROR problem\n";
#endif
				continue;
			}
			vid[i].push_back(it->first);
			val[i].push_back(it->second);
		}
	}
	if (vid.size() > 0)
		solver.append_orientation_cons(vid, val);

}

void Optimizer::append_orientation_cons_curve_based(map<std::pair<int, int>, int>& mp, Mosek_c_wrapper &solver)
{
	//ccid,edgeid,confid-x
	std::vector<std::vector<std::vector<int>>> ori_plus, ori_minus;
	for (int i = 0; i < ori_cons_with_confid.size(); i++) {
		std::vector<std::vector<int>> ori_plus_ci, ori_minus_ci;
		std::map<int, map<int, std::vector<orientation_Expr_with_id>>>::iterator it = ori_cons_with_confid[i].begin();
		for (; it != ori_cons_with_confid[i].end(); it++) {//go through all edges
			int curveid = it->first;

			map<int, std::vector<orientation_Expr_with_id>>::iterator diff_incident_faces_it = it->second.begin();
			for (; diff_incident_faces_it != it->second.end(); diff_incident_faces_it++) {//each edge may have differnt incident faces as ci

				std::vector<int> ori_plus_ci_edgeit, ori_minus_ci_edgeit;
				bool pushed1 = false, pushed2 = false;

				//std::map<int, int>::iterator nit = node.find(edgeid);

				//not fixed
				for (int j = 0; j < diff_incident_faces_it->second.size(); j++) {
					int confid = (diff_incident_faces_it->second)[j].confid;
					bool dir = (diff_incident_faces_it->second)[j].positive_flg;


					if (dir) {
						//positive
						ori_plus_ci_edgeit.push_back(mp[std::make_pair(curveid, confid)]);
						pushed1 = true;

					}
					else {
						//negative
						ori_minus_ci_edgeit.push_back(mp[std::make_pair(curveid, confid)]);
						pushed2 = true;
					}
				}



				if (pushed1)
					ori_plus_ci.push_back(ori_plus_ci_edgeit);
				if (pushed2)
					ori_minus_ci.push_back(ori_minus_ci_edgeit);
			}


		}

		ori_plus.push_back(ori_plus_ci);
		ori_minus.push_back(ori_minus_ci);
	}


	//orientation constraints
	//if dual edge induced direction follow component normal, then the contraint should be pushed to positive slot; if not, negative
	//the dual edge positive rely on the direction of incident face with smaller incident idx in the incident list
	int cons_counter = 0;
	std::vector<map<int, double>> trilist;//use trilist to combine duplicate term.//rowid-colid-val
	for (int m = 0; m < ori_plus.size(); m++) {
		int siz = ori_plus[m].size();
		for (int i = 0; i < siz; i++) {
			int j = (i + 1) % siz;
			trilist.push_back(map<int, double>{});
			int consid = trilist.size();
			consid--;


			for (int ss = 0; ss < ori_plus[m][i].size(); ss++) {
				auto ret = trilist[consid].insert(std::make_pair(ori_plus[m][i][ss], 1.0));
				if (ret.second == false) {
					//already exist
					ret.first->second += 1.0;
				}

			}
			for (int ss = 0; ss < ori_plus[m][j].size(); ss++) {
				auto ret = trilist[consid].insert(std::make_pair(ori_plus[m][j][ss], -1.0));
				if (ret.second == false) {
					//already exist
					ret.first->second += -1.0;
				}
			}

		}
	}

	for (int m = 0; m < ori_minus.size(); m++) {
		int siz = ori_minus[m].size();
		for (int i = 0; i < siz; i++) {
			int j = (i + 1) % siz;
			trilist.push_back(map<int, double>{});
			int consid = trilist.size();
			consid--;



			for (int ss = 0; ss < ori_minus[m][i].size(); ss++) {
				auto ret = trilist[consid].insert(std::make_pair(ori_minus[m][i][ss], 1.0));
				if (ret.second == false) {
					//already exist
					ret.first->second += 1.0;
				}

			}
			for (int ss = 0; ss < ori_minus[m][j].size(); ss++) {
				auto ret = trilist[consid].insert(std::make_pair(ori_minus[m][j][ss], -1.0));
				if (ret.second == false) {
					//already exist
					ret.first->second += -1.0;
				}
			}

		}
	}

	std::vector<std::vector<int>> vid;
	std::vector<std::vector<double>> val;

	for (int i = 0; i < trilist.size(); i++) {
		vid.push_back(std::vector<int>{});
		val.push_back(std::vector<double>{});
		for (auto it = trilist[i].begin(); it != trilist[i].end(); it++) {
			if (abs(it->second - 0) < 1e-7) {
#ifndef _DISABLE_FILLING_WARNING
				cout << "problem,";
#endif
#if defined(_RECORD_ERROR_LOGS)
				//(*_log_off) << _file_under_processing << " ERROR problem\n";
#endif
				continue;
			}
			vid[i].push_back(it->first);
			val[i].push_back(it->second);
		}
	}
	if (vid.size() > 0)
		solver.append_orientation_cons(vid, val);

}


void Optimizer::append_x_regular_cons(map<std::pair<int, int>, int>& mp, Mosek_c_wrapper &solver)
{
	std::vector<std::vector<int>> vid;
	std::vector<std::vector<double>> val;

	for (int i = 0; i < _gh.edge_list.size(); i++) {
		int degree = _gh.edge_list[i].incident_components.size();

		vid.push_back(std::vector<int>{});
		val.push_back(std::vector<double>{});

		for (int j = 0; j < configuration_list_all_degrees[degree - 2].size(); j++) {
			auto ret = mp.find(std::make_pair(i, j));
			if (ret != mp.end()) {
				vid[i].push_back(ret->second);
				val[i].push_back(1.0);
			}
		}


	}
	if (vid.size() > 0)
		solver.append_regular_cons(vid, val);
}
void Optimizer::append_x_regular_cons_curve_based(map<std::pair<int, int>, int>& mp, Mosek_c_wrapper &solver)
{
	std::vector<std::vector<int>> vid;
	std::vector<std::vector<double>> val;

	for (int i = 0; i < _gh.non_manifold_curves.size(); i++) {
		int degree = _gh.edge_list[_gh.non_manifold_curves[i][0]].incident_components.size();

		vid.push_back(std::vector<int>{});
		val.push_back(std::vector<double>{});

		for (int j = 0; j < configuration_list_all_degrees[degree - 2].size(); j++) {
			auto ret = mp.find(std::make_pair(i, j));
			if (ret != mp.end()) {
				vid[i].push_back(ret->second);
				val[i].push_back(1.0);
			}
		}


	}
	if (vid.size() > 0)
		solver.append_regular_cons(vid, val);
}

void Optimizer::append_global_cons_with_all_fixed(map<std::pair<int, int>, int> &mp, std::vector<double> solution, Mosek_c_wrapper *solver)//solution contain do not s
{
	//discard
}

void Optimizer::append_fixed_variables(map<std::pair<int, int>, int>& mp, Mosek_c_wrapper & solver, std::map<int, int> node)
{// as \sum x=1, so only push the varible=1 cons
	std::vector<std::vector<int>> vid;
	std::vector<std::vector<double>> val;

	for (auto it = node.begin(); it != node.end(); it++) {
		int edgeid = it->first;
		int confid = it->second;
		vid.push_back(std::vector<int>{mp[std::make_pair(edgeid, confid)]});
		val.push_back(std::vector<double>{1.0});
	}

	if (vid.size() > 0)
		solver.append_fixed_var_cons(vid, val);
}

void Optimizer::append_fixed_variables_debug_use(Mosek_c_wrapper & solver, std::vector<double> &s)
{
	std::vector<std::vector<int>> vid;
	std::vector<std::vector<double>> val;

	for (int i = 0; i < s.size(); i++) {

		vid.push_back(std::vector<int>{i + 1});
		val.push_back(std::vector<double>{1.0});
	}

	if (vid.size() > 0)
		solver.append_fixed_var_cons_debug_use(vid, val, s);
}


std::array<double, 5> Optimizer::solve_problem_with_fixed_solution_capi(std::map<int, int>& node, std::vector<std::vector<double>>& result_solution)
{
	result_solution.clear();

	Mosek_c_wrapper msk_solver;;

	msk_solver.cloned_from(base_msk_solver);

	append_fixed_variables(var_map, msk_solver, node);

	msk_solver.solve();

	double s_bound = 0;
	std::vector<double> solution;
	msk_solver.get_solution(s_bound, solution, n_opt_vars);


	for (int i = 0; i < _gh.non_manifold_curves.size(); i++) {
		int degree = _gh.edge_list[_gh.non_manifold_curves[i][0]].incident_components.size();
		result_solution.push_back(std::vector<double>{});
		for (int j = 0; j < configuration_list_all_degrees[degree - 2].size(); j++) {
			result_solution[i].push_back(solution[var_map[std::make_pair(i, j)] - 1]);
		}
	}

	double visual_sc = 0;
	msk_solver.~Mosek_c_wrapper();

	for (int i = 0; i < obj_vid.size(); i++) {
		visual_sc += result_solution[obj_vid[i].first][obj_vid[i].second] * obj_val[i];
	}


#ifdef _CHECK_EACH_ENERGY_TERM
	double unary_s = 0.0;
	for (auto it = valid_val_unary.begin(); it != valid_val_unary.end(); it++) {
		int curvid = it->first.first;
		int confid = it->first.second;
		double val = it->second;
		unary_s += result_solution[curvid][confid] * val;
	}

	double binary_s = 0.0;
	for (auto it = valid_val_binary.begin(); it != valid_val_binary.end(); it++) {
		int curvid = it->first.first;
		int confid = it->first.second;
		double val = it->second;
		binary_s += result_solution[curvid][confid] * val;
	}


#ifdef _OUTPUT_DEBUG_LOGS
	std::vector<double> unary_term(_gh.node_list.size() * 2, 0.0);
	for (int i = 0; i < unary_score.size(); i++) {
		//accumulate prob of conf
		{
			double xprob = 0.0;
			for (auto it = unary_score[i].confid_pos.begin();
				it != unary_score[i].confid_pos.end(); it++) {

				xprob += result_solution[it->first][it->second];
			}

			unary_term[2 * i] = unary_score[i].pos_score*xprob*unary_score[i].weight;

		}

		{
			double xprob = 0.0;
			for (auto it = unary_score[i].confid_neg.begin();
				it != unary_score[i].confid_neg.end(); it++) {

				xprob += result_solution[it->first][it->second];
			}

			unary_term[2 * i + 1] = unary_score[i].neg_score*xprob*unary_score[i].weight;

		}
	}



	//unary from binary
	for (auto it = unary_score_from_binary.begin();
		it != unary_score_from_binary.end(); it++) {
		int curveid = it->first.first;
		int dualedgeid = it->first.second;

		int someedgeid = _gh.non_manifold_curves[curveid][0];
		int ccidi = _gh.edge_list[someedgeid].incident_components[it->second.inidi];
		int ccidj = _gh.edge_list[someedgeid].incident_components[it->second.inidj];


		if (ccidi != ccidj) {
			std::cout << "Error\n";
#if defined(_RECORD_ERROR_LOGS)
			(*_log_off) << _file_under_processing << " ERROR in solving\n";
			(*_log_off).flush();
#endif
		}

		double xprob = 0;
		for (auto ti = it->second.confid.begin();
			ti != it->second.confid.end(); ti++) {
			int confid = *ti;
			xprob += result_solution[curveid][confid];
		}


		double sc = it->second.weight*it->second.score*xprob;


		unary_term[2 * ccidi + it->second.diri] += sc;
	}
#endif

	return{ s_bound*_sWeight + visual_sc,s_bound*_sWeight,visual_sc, unary_s, binary_s };
#else
	return{ s_bound + visual_sc,s_bound,visual_sc, 0.0, 0.0 };
#endif


}

//rounding
void  Optimizer::maximal_rounding(std::vector<std::vector<double>>& input, std::vector<std::vector<double>>& output) {
	output.clear();

	for (int i = 0; i < input.size(); i++) {
		double maxval = -1;
		int maxid = -1;

		for (int j = 0; j < input[i].size(); j++) {
			if (maxval < input[i][j]) {
				maxval = input[i][j];
				maxid = j;
			}
		}

		std::vector<double> rounded_vec;
		rounded_vec.resize(input[i].size(), 0);
		rounded_vec[maxid] = 1;

		output.push_back(rounded_vec);
	}
}


//feassible

bool Optimizer::if_x_feasible(std::vector<std::vector<double>>& x) {

	bool res = if_sum_x_equal_to_one(x);

	if (res) {
		res = if_x_meet_orientation_cons(x);
	}

	return res;
}

#define _FEASIABLE_EPS_ 1e-5
bool Optimizer::if_sum_x_equal_to_one(std::vector<std::vector<double>>& x) {
	for (int i = 0; i < x.size(); i++) {
		double sum_res = 0;
		for (int j = 0; j < x[i].size(); j++) {
			sum_res += x[i][j];
		}

		if (abs(sum_res - 1.0) > _FEASIABLE_EPS_) return false;
	}
	return true;
}

bool Optimizer::if_branched_solution_meet_cons(std::map<int, int> node) {
	std::vector<std::vector<double>> ori_plus, ori_minus;
	for (int i = 0; i < ori_cons_with_confid.size(); i++) {
		std::vector<double> ori_plus_ci, ori_minus_ci;
		std::map<int, map<int, std::vector<orientation_Expr_with_id>>>::iterator it = ori_cons_with_confid[i].begin();
		for (; it != ori_cons_with_confid[i].end(); it++) {//go through all edges
			int edgeid = it->first;

			map<int, std::vector<orientation_Expr_with_id>>::iterator diff_incident_faces_it = it->second.begin();
			for (; diff_incident_faces_it != it->second.end(); diff_incident_faces_it++) {//each edge may have differnt incident faces as ci

				double ori_plus_ci_edgeit, ori_minus_ci_edgeit;
				bool pushed1 = false, pushed2 = false;

				int picked_confid = -1;
				if (node.find(edgeid) == node.end()) {
					continue;
				}
				else {
					picked_confid = node[edgeid];
				}

				for (int j = 0; j < diff_incident_faces_it->second.size(); j++) {
					int confid = (diff_incident_faces_it->second)[j].confid;
					bool dir = (diff_incident_faces_it->second)[j].positive_flg;
					if (confid == picked_confid) {
						pushed1 = true;
						pushed2 = true;
						if (dir) {
							ori_plus_ci_edgeit = 1.0;
							ori_minus_ci_edgeit = 0.0;
						}
						else {
							ori_plus_ci_edgeit = 0.0;
							ori_minus_ci_edgeit = 1.0;
						}
					}
				}
				if (pushed1)
					ori_plus_ci.push_back(ori_plus_ci_edgeit);
				if (pushed2)
					ori_minus_ci.push_back(ori_minus_ci_edgeit);
			}


		}

		ori_plus.push_back(ori_plus_ci);
		ori_minus.push_back(ori_minus_ci);
	}

	for (int m = 0; m < ori_plus.size(); m++) {
		int siz = ori_plus[m].size();
		for (int i = 0; i < siz; i++) {
			int j = (i + 1) % siz;
			if (abs(ori_plus[m][i] - ori_plus[m][j]) > 1e-7)return false;
			if (abs(ori_minus[m][i] - ori_minus[m][j]) > 1e-7)return false;

		}
	}

	return true;
}

bool Optimizer::if_x_meet_orientation_cons(std::vector<std::vector<double>>& x) {
	std::vector<std::vector<double>> ori_plus, ori_minus;
	
	for (int i = 0; i < ori_cons_with_confid.size(); i++) {
		std::vector<double> ori_plus_ci, ori_minus_ci;
		std::map<int, map<int, std::vector<orientation_Expr_with_id>>>::iterator it = ori_cons_with_confid[i].begin();
		for (; it != ori_cons_with_confid[i].end(); it++) {//go through all edges
			int edgeid = it->first;

			map<int, std::vector<orientation_Expr_with_id>>::iterator diff_incident_faces_it = it->second.begin();
			for (; diff_incident_faces_it != it->second.end(); diff_incident_faces_it++) {//each edge may have differnt incident faces as ci

				double ori_plus_ci_edgeit, ori_minus_ci_edgeit;
				bool pushed1 = false, pushed2 = false;

				int picked_confid = -1;
				for (int ss = 0; ss < x[edgeid].size(); ss++) {
					if (abs(x[edgeid][ss] - 1.0) < 1e-7) {
						picked_confid = ss;
						break;
					}
				}

				for (int j = 0; j < diff_incident_faces_it->second.size(); j++) {
					int confid = (diff_incident_faces_it->second)[j].confid;
					bool dir = (diff_incident_faces_it->second)[j].positive_flg;
					if (confid == picked_confid) {
						pushed1 = true;
						pushed2 = true;
						if (dir) {
							ori_plus_ci_edgeit = 1.0;
							ori_minus_ci_edgeit = 0.0;
						}
						else {
							ori_plus_ci_edgeit = 0.0;
							ori_minus_ci_edgeit = 1.0;
						}
					}
				}
				if (pushed1)
					ori_plus_ci.push_back(ori_plus_ci_edgeit);
				if (pushed2)
					ori_minus_ci.push_back(ori_minus_ci_edgeit);
			}


		}

		ori_plus.push_back(ori_plus_ci);
		ori_minus.push_back(ori_minus_ci);
	}



	for (int m = 0; m < ori_plus.size(); m++) {
		int siz = ori_plus[m].size();
		for (int i = 0; i < siz; i++) {
			int j = (i + 1) % siz;
			if (abs(ori_plus[m][i] - ori_plus[m][j]) > 1e-7)return false;
			if (abs(ori_minus[m][i] - ori_minus[m][j]) > 1e-7)return false;

		}
	}

	return true;
}

//energy
double Optimizer::algebraic_conn_bound(std::vector<std::vector<double>>& x) {

	return 0.0;

}

std::array<double, 5> Optimizer::algebraic_conn_bound_capi(std::vector<std::vector<double>>& x)
{

	double s_bound = -1.0;

	double vis_sc = 0.0;

	for (int i = 0; i < obj_vid.size(); i++) {
		vis_sc += x[obj_vid[i].first][obj_vid[i].second] * obj_val[i];
	}


#ifdef _CHECK_EACH_ENERGY_TERM
	double unary_s = 0.0;
	for (auto it = valid_val_unary.begin(); it != valid_val_unary.end(); it++) {
		int curvid = it->first.first;
		int confid = it->first.second;
		double val = it->second;
		unary_s += x[curvid][confid] * val;
	}

	double binary_s = 0.0;
	for (auto it = valid_val_binary.begin(); it != valid_val_binary.end(); it++) {
		int curvid = it->first.first;
		int confid = it->first.second;
		double val = it->second;
		binary_s += x[curvid][confid] * val;
	}
	return{ vis_sc + s_bound*_sWeight, s_bound*_sWeight, vis_sc,unary_s ,binary_s };
#else
	return{ vis_sc + s_bound, s_bound, vis_sc,0,0 };
#endif

}

std::array<double, 5> Optimizer::algebraic_conn_bound_capi_debug_use(std::vector<std::vector<double>>& x)
{
	Mosek_c_wrapper msk_solver;;

	msk_solver.cloned_from(base_msk_solver);

	std::vector<double> s;
	for (int i = 0; i < x.size(); i++) {
		for (int j = 0; j < x[i].size(); j++) {
			s.push_back(x[i][j]);
		}
	}

	append_fixed_variables_debug_use(msk_solver, s);

	bool succ = msk_solver.solve();
	if (!succ) {
		return{ -1,-1,-1,-1,-1 };
	}

	double s_bound = 0;
	std::vector<double> solution;
	msk_solver.get_solution(s_bound, solution, n_opt_vars);

	msk_solver.~Mosek_c_wrapper();

	double vis_sc = 0.0;

	for (int i = 0; i < obj_vid.size(); i++) {
		vis_sc += x[obj_vid[i].first][obj_vid[i].second] * obj_val[i];
	}


#ifdef _CHECK_EACH_ENERGY_TERM
	double unary_s = 0.0;
	for (auto it = valid_val_unary.begin(); it != valid_val_unary.end(); it++) {
		int curvid = it->first.first;
		int confid = it->first.second;
		double val = it->second;
		unary_s += x[curvid][confid] * val;
	}

	double binary_s = 0.0;
	for (auto it = valid_val_binary.begin(); it != valid_val_binary.end(); it++) {
		int curvid = it->first.first;
		int confid = it->first.second;
		double val = it->second;
		binary_s += x[curvid][confid] * val;
	}
	return{ vis_sc + s_bound*_sWeight, s_bound*_sWeight, vis_sc,unary_s ,binary_s };
#else
	return{ vis_sc + s_bound, s_bound, vis_sc,0,0 };
#endif
}

void Optimizer::update_unary_and_binary_score(double CUTTING_EPS)
{
	Tribunal tri;
	_msh->update_mesh_properties();
	tri.setMesh(_msh);
	tri.setPrimalDualGraph(&_gh);
	tri.setSamplig(&_sampler);
	tri.setScoreList(&(binary_score), &(unary_score_from_binary), &(unary_score));
	tri.init();
	tri.update_unary_score(CUTTING_EPS);
	tri.upate_binary_score();
	tri.update_weight(non_manifold_curve_length, _uWeight, _bWeight);

}

void Optimizer::append_objective_with_visual_score(map<std::pair<int, int>, int>& mp, Mosek_c_wrapper & solver)
{
#ifdef _CHECK_EACH_ENERGY_TERM
	valid_val_unary.clear();
	valid_val_binary.clear();
#endif

	//maximize!!
	std::map<std::pair<int, int>, double> valid_val;
	//unary term
	for (int i = 0; i < unary_score.size(); i++) {
		{//pos
			double score = unary_score[i].pos_score;
			auto confs = unary_score[i].confid_pos;
			double weight = unary_score[i].weight;

			for (auto it = confs.begin(); it != confs.end(); it++) {
				int curveid = it->first;
				int confid = it->second;
				

				//auto nid = mp.find(std::make_pair(edgeid, confid));
				auto fit = valid_val.insert(std::make_pair(std::make_pair(curveid, confid), 0.0));
				fit.first->second += score*weight;// add score to the term
#ifdef _CHECK_EACH_ENERGY_TERM
				auto fit2 = valid_val_unary.insert(std::make_pair(std::make_pair(curveid, confid), 0.0));
				fit2.first->second += score*weight;// add score to the term
#endif
			}
		}
		{//neg
			double score = unary_score[i].neg_score;
			auto confs = unary_score[i].confid_neg;
			double weight = unary_score[i].weight;

			for (auto it = confs.begin(); it != confs.end(); it++) {
				int curveid = it->first;
				int confid = it->second;


				//auto nid = mp.find(std::make_pair(edgeid, confid));
				auto fit = valid_val.insert(std::make_pair(std::make_pair(curveid, confid), 0.0));
				fit.first->second += score*weight;// add score to the term
#ifdef _CHECK_EACH_ENERGY_TERM
				auto fit2 = valid_val_unary.insert(std::make_pair(std::make_pair(curveid, confid), 0.0));
				fit2.first->second += score*weight;// add score to the term
#endif
			}
		}
	}


	//unary term from binary
	for (auto it = unary_score_from_binary.begin(); it != unary_score_from_binary.end(); it++) {
		int curveid = it->first.first;
		double score = it->second.score;
		double weight = it->second.weight;

		for (auto iit = it->second.confid.begin(); iit != it->second.confid.end(); iit++) {
			int confid = *iit;
			//auto nid = mp.find(std::make_pair(curveid, confid));
			auto fit = valid_val.insert(std::make_pair(std::make_pair(curveid, confid), 0.0));
			fit.first->second += score*weight;// add score to the term
#ifdef _CHECK_EACH_ENERGY_TERM
			auto fit2 = valid_val_unary.insert(std::make_pair(std::make_pair(curveid, confid), 0.0));
			fit2.first->second += score*weight;// add score to the term
#endif
		}
	}

	//binary term
	for (auto it = binary_score.begin(); it != binary_score.end(); it++) {
		int curveid = it->first.first;
		double score = it->second.score;
		double weight = it->second.weight;
		for (auto iit = it->second.confid.begin(); iit != it->second.confid.end(); iit++) {
			int confid = *iit;
			//auto nid = mp.find(std::make_pair(curveid, confid));
			auto fit = valid_val.insert(std::make_pair(std::make_pair(curveid, confid), 0.0));
			fit.first->second += score*weight;// add score to the term	
#ifdef		_CHECK_EACH_ENERGY_TERM
			auto fit2 = valid_val_binary.insert(std::make_pair(std::make_pair(curveid, confid), 0.0));
			fit2.first->second += score*weight;// add score to the term
#endif
		}
	}

	//turn into vector
	std::vector<int> vid;
	//std::vector<double> val;
	obj_vid.clear();
	obj_val.clear();
	for (auto it = valid_val.begin(); it != valid_val.end(); it++) {
		if (abs(it->second) > 1e-7) {
			auto iit = mp.find(it->first);
			vid.push_back(iit->second);
			obj_vid.push_back(it->first);
			obj_val.push_back(it->second);
		}
	}

	if (vid.size() > 0)
		solver.append_objective_with_visual_score(vid, obj_val, _sWeight);
}


//write mesh
void Optimizer::get_directioend_component_list(std::vector<int> &x, std::vector<int> &ccdir, std::vector<std::vector<int>> &ncclist) {
	ccdir.clear();
	ccdir.resize(_gh.node_list.size(), -1);
	ncclist.clear();

	std::vector<bool> if_component_visited;
	if_component_visited.resize(_gh.node_list.size(), false);

	for (int i = 0; i < if_component_visited.size(); i++) {


		if (if_component_visited[i])continue;

		//not visited
		std::vector<int> cc_list;
		std::vector<int> queue;//ccid, true-postivie

		queue.push_back(i);
		if_component_visited[i] = true;

		while (queue.size() != 0) {
			int qsiz = queue.size();
			int current_ccid = queue[qsiz - 1];
			//bool current_dir= queue[qsiz - 1].second;
			queue.pop_back();


			cc_list.push_back(current_ccid);

			if (_gh.node_list[current_ccid].incident_edges.size() == 0) {
				ccdir[current_ccid]=0;
				continue;
			}

			for (int j = 0; j < _gh.node_list[current_ccid].incident_edges.size(); j++) {

				//the jth edge of component i
				int edgeid = _gh.node_list[current_ccid].incident_edges[j];

				int picked_conf = x[edgeid];
				//check if the conf really connect c
				int degree = _gh.edge_list[edgeid].incident_components.size();
				std::vector<int> picked_dual_edges = configuration_list_all_degrees[degree - 2][picked_conf];
				//check each dual edge 
				for (int m = 0; m < picked_dual_edges.size(); m++) {
					int idi = dual_edge_list_all_degrees[degree - 2][picked_dual_edges[m]].fi;
					int idj = dual_edge_list_all_degrees[degree - 2][picked_dual_edges[m]].fj;
					int dir = dual_edge_list_all_degrees[degree - 2][picked_dual_edges[m]].dir;
					//int dir = dual_edge_list_all_degrees[degree - 2][picked_dual_edges[m]].dir;
					if (idi == -1 || idj == -1) {
						int idd = idi == -1 ? idj : idi;
						int ccid = _gh.edge_list[edgeid].incident_components[idd];

#ifdef _OUTPUT_DEBUG_LOGS
						if (ccdir[ccid] != -1 && ccdir[ccid] != dir) {
							std::cout << "Error";
#if defined(_RECORD_ERROR_LOGS)
							(*_log_off) << _file_under_processing << " ERROR get_directioend_component_list\n";
							(*_log_off).flush();
#endif
						}
#endif
						ccdir[ccid] = dir;

						continue;
					}

					int fidi = _gh.edge_list[edgeid].incident_components[idi];
					int fidj = _gh.edge_list[edgeid].incident_components[idj];


					//original face id
					int fidfi = _msh->GetEdges()[_gh.edge_list[edgeid].original_edge_id].incident_faces[idi];
					int fidfj = _msh->GetEdges()[_gh.edge_list[edgeid].original_edge_id].incident_faces[idj];

					BlackMesh::BlackMesh<double>::Triangle tr1 = _msh->GetTriangles()[fidfi];
					BlackMesh::BlackMesh<double>::Triangle tr2 = _msh->GetTriangles()[fidfj];

					int flg = if_face_coherant_based_on_edge(_msh->GetEdges()[_gh.edge_list[edgeid].original_edge_id].vertices.first,
						_msh->GetEdges()[_gh.edge_list[edgeid].original_edge_id].vertices.second,
						tr1, tr2);

					//current face flipped, coherant, conencted flipped, p-keep original dir, n-flip
					//p,p,p,
					//p,n,n,
					//n,p,n,
					//n,n,p,

					//the dual edge do connect ci
					if (fidi == current_ccid) {
						//change the dir of fidi and fidj
#ifdef _OUTPUT_DEBUG_LOGS
						if (ccdir[fidi] != -1 && ccdir[fidi] != dir) {
							std::cout << "Error";
#if defined(_RECORD_ERROR_LOGS)
							(*_log_off) << _file_under_processing << " ERROR get_directioend_component_list\n";
							(*_log_off).flush();
#endif
						}
#endif
						ccdir[fidi] = dir;
						if (flg == 1)
						{
#ifdef _OUTPUT_DEBUG_LOGS
							if (ccdir[fidj] != -1 && ccdir[fidj] != dir) {
								std::cout << "Error";
#if defined(_RECORD_ERROR_LOGS)
								(*_log_off) << _file_under_processing << " ERROR get_directioend_component_list\n";
								(*_log_off).flush();
#endif
							}
#endif
							ccdir[fidj] = dir;
						}
						else {
#ifdef _OUTPUT_DEBUG_LOGS
							if (ccdir[fidj] != -1 && ccdir[fidj] != (dir + 1) % 2) {
								std::cout << "Error";
#if defined(_RECORD_ERROR_LOGS)
								(*_log_off) << _file_under_processing << " ERROR get_directioend_component_list\n";
								(*_log_off).flush();
#endif
							}
#endif
							ccdir[fidj] = (dir + 1) % 2;
						}

						//should push fidj to queue
						if (!if_component_visited[fidj]) {
							if_component_visited[fidj] = true;
							queue.push_back(fidj);
						}
					}
					if (fidj == current_ccid) {
						//change the dir of ccidi and ccidj
#ifdef _OUTPUT_DEBUG_LOGS
						if (ccdir[fidi] != -1 && ccdir[fidi] != dir) {
							std::cout << "Error";
#if defined(_RECORD_ERROR_LOGS)
							(*_log_off) << _file_under_processing << " ERROR get_directioend_component_list\n";
							(*_log_off).flush();
#endif
						}
#endif
						ccdir[fidi] = dir;
						if (flg == 1)
						{
#ifdef _OUTPUT_DEBUG_LOGS
							if (ccdir[fidj] != -1 && ccdir[fidj] != dir) {
								std::cout << "Error";
#if defined(_RECORD_ERROR_LOGS)
								(*_log_off) << _file_under_processing << " ERROR get_directioend_component_list\n";
								(*_log_off).flush();
#endif
							}
#endif
							ccdir[fidj] = dir;
						}
						else {
#ifdef _OUTPUT_DEBUG_LOGS
							if (ccdir[fidj] != -1 && ccdir[fidj] != (dir + 1) % 2) {
								std::cout << "Error";
#if defined(_RECORD_ERROR_LOGS)
								(*_log_off) << _file_under_processing << " ERROR get_directioend_component_list\n";
								(*_log_off).flush();
#endif
							}
#endif
							ccdir[fidj] = (dir + 1) % 2;
						}

						//should push fidi to queue
						if (!if_component_visited[fidi]) {
							if_component_visited[fidi] = true;
							queue.push_back(fidi);
						}
					}
				}
			}
		}
		ncclist.push_back(cc_list);
	}
}

void Optimizer::change_mesh_topo_with_solution(std::vector<int>& x, std::vector<bool> & del_cc_ls)
{
	//x is edge-based
	auto build_oriedge_to_graphedge_map = [](map<int, int> &res, BlackMesh::BlackMesh<double> *msh, Primal_Dual_graph *gh) {
		res.clear();

		for (int i = 0; i < gh->edge_list.size(); i++) {

			res.insert(std::make_pair(gh->edge_list[i].original_edge_id, i));

		}

#ifdef _OUTPUT_DEBUG_LOGS
		for (int i = 0; i < gh->edge_list.size(); i++) {
			if (res[gh->edge_list[i].original_edge_id] != i) {
				std::cout << "The map between edges is not matched\n";
#if defined(_RECORD_ERROR_LOGS)
				(*_log_off) << _file_under_processing << " ERROR write_mesh_change_topo_with_solution\n";
				(*_log_off).flush();
#endif
			}
		}
#endif
	};

	auto get_incident_non_manifold_edges = [](map<int, int> &mp, std::vector<int> &ori_edgeid, std::vector<int> &gh_edgeid) {
		gh_edgeid.clear();
		for (int i = 0; i < ori_edgeid.size(); i++) {
			auto it = mp.find(ori_edgeid[i]);
			if (it == mp.end()) {
				//do nothing
			}
			else {
				gh_edgeid.push_back(it->second);
			}
		}
	};

	auto get_eq_face_set_on_vtx = [&](int vertexid, const map<int, int> &msh2ghmap, const std::vector<int>& x, std::vector<std::vector<int>> &output_groups) {
		std::set<int> visited_face;
		std::set<int> all_face;
		//get all faces incident to vtx
		auto incident_edges_to_vtx = _msh->GetVertices()[vertexid].incident_edges;
		for (int i = 0; i < incident_edges_to_vtx.size(); i++) {
			all_face.insert(_msh->GetEdges()[incident_edges_to_vtx[i]].incident_faces.begin(),
				_msh->GetEdges()[incident_edges_to_vtx[i]].incident_faces.end());
		}

		//go through all faces
		output_groups.clear();
		for (auto it = all_face.begin(); it != all_face.end(); it++) {
			int faceid = *it;
			if (visited_face.find(faceid) != visited_face.end())continue;

			std::queue<int> list;
			list.push(faceid);
			visited_face.insert(faceid);
			output_groups.push_back(std::vector<int>{faceid});
			int siz = output_groups.size();
			int groupid = siz - 1;

			while (list.size() != 0) {
				int current_face = list.front();
				list.pop();

				for (int i = 0; i < _msh->GetTriangles()[current_face].incident_edges.size(); i++) {
					int edgeidmsh = _msh->GetTriangles()[current_face].incident_edges[i];

					if (_msh->GetEdges()[edgeidmsh].is_nonmanifold_edge == false) {
						//not nedge
						for (int j = 0; j < _msh->GetEdges()[edgeidmsh].incident_faces.size(); j++) {
							int query_face_id = _msh->GetEdges()[edgeidmsh].incident_faces[j];
							if (visited_face.find(query_face_id) != visited_face.end())continue;
							//if query face have vtx
							bool hit = false;

							if (all_face.find(query_face_id) != all_face.end())hit = true;

							if (hit) {
								list.push(query_face_id);
								visited_face.insert(query_face_id);
								output_groups[groupid].push_back(query_face_id);
							}
						}

					}
					else {
						//nedge
						int edgeidgh = msh2ghmap.find(edgeidmsh)->second;
						int degree = _gh.edge_list[edgeidgh].incident_components.size();
						auto pickedconf = configuration_list_all_degrees[degree - 2][x[edgeidgh]];

						//get the incident id of current face
						int incident_id_of_current_face = -1;
						for (int j = 0; j < _msh->GetEdges()[edgeidmsh].incident_faces.size(); j++) {
							if (_msh->GetEdges()[edgeidmsh].incident_faces[j] == current_face) {
								incident_id_of_current_face = j;
								break;
							}
						}

						//find the dual edge picked that have such incident id
						int anotherId = -2;
						for (int j = 0; j < pickedconf.size(); j++) {
							auto picked_dualedge = dual_edge_list_all_degrees[degree - 2][pickedconf[j]];

							if (picked_dualedge.fi == incident_id_of_current_face) {
								if (anotherId != -2) {
									std::cout << "ERROR HERR\n";
#if defined(_RECORD_ERROR_LOGS)
									(*_log_off) << _file_under_processing << " ERROR write_mesh_change_topo_with_solution\n";
									(*_log_off).flush();
#endif
								}
								anotherId = picked_dualedge.fj;
							}
							else if (picked_dualedge.fj == incident_id_of_current_face) {
								if (anotherId != -2) {
									std::cout << "ERROR HERR\n";
#if defined(_RECORD_ERROR_LOGS)
									(*_log_off) << _file_under_processing << " ERROR write_mesh_change_topo_with_solution\n";
									(*_log_off).flush();
#endif
								}
								anotherId = picked_dualedge.fi;
							}
						}

						if (anotherId == -2) {
							std::cout << "ERROR HERR\n";
#if defined(_RECORD_ERROR_LOGS)
							(*_log_off) << _file_under_processing << " ERROR write_mesh_change_topo_with_solution\n";
							(*_log_off).flush();
#endif
						}
						else if (anotherId == -1) {
							//current face is break to the nedge
						}
						else {
							int query_face_id = _msh->GetEdges()[edgeidmsh].incident_faces[anotherId];
							if (visited_face.find(query_face_id) != visited_face.end())continue;
							//if query face have vtx
							bool hit = false;
							if (all_face.find(query_face_id) != all_face.end())hit = true;
							if (hit) {
								list.push(query_face_id);
								visited_face.insert(query_face_id);
								output_groups[groupid].push_back(query_face_id);
							}
						}

					}

				}

			}
		}

	};

	std::map<int, int> old_to_new_edge;
	build_oriedge_to_graphedge_map(old_to_new_edge, _msh, &_gh);

	std::map<int, std::vector<std::vector<int>>> vertex_split_groups;
	for (int i = 0; i < _gh.edge_list.size(); i++) {

		int vtx;


		vtx = _msh->GetEdges()[_gh.edge_list[i].original_edge_id].vertices.first;
		auto it = vertex_split_groups.insert(std::make_pair(vtx, std::vector<std::vector<int>>{}));
		if (it.second == false) {
			//exist
		}
		else {
			get_eq_face_set_on_vtx(vtx, old_to_new_edge, x, it.first->second);
		}



		vtx = _msh->GetEdges()[_gh.edge_list[i].original_edge_id].vertices.second;
		it = vertex_split_groups.insert(std::make_pair(vtx, std::vector<std::vector<int>>{}));
		if (it.second == false) {
			//exist
		}
		else {
			get_eq_face_set_on_vtx(vtx, old_to_new_edge, x, it.first->second);
		}
	}

	/*****/
	//get component list
	std::vector<int> ccdir;
	std::vector<std::vector<int>> ncclist;
	get_directioend_component_list(x, ccdir, ncclist);

	BlackMesh::BlackMesh<double> *nmesh = new BlackMesh::BlackMesh<double>();

	for (int i = 0; i < _msh->GetNumVertices(); i++) {
		nmesh->insert_vtx(_msh->GetVerticesEditable()[i].pos);
	}

	//write extra vtx
	std::map<std::pair<int, int>, int> idx_mapping;//ori_vtxid,ccid-newvtxid
	int idcounter = _msh->GetNumVertices();
	idcounter--;
	for (auto it = vertex_split_groups.begin(); it != vertex_split_groups.end(); it++) {
		//0th, push original id

		int oriid = it->first;
		int groupid = 0;
		idx_mapping.insert(std::make_pair(std::make_pair(oriid, 0), oriid));


		for (int i = 1; i < it->second.size(); i++) {
			int oriid = it->first;
			//off << "v " << _msh->GetVertices()[oriid].pos[0]
			//	<< " " << _msh->GetVertices()[oriid].pos[1]
			//	<< " " << _msh->GetVertices()[oriid].pos[2] << "\n";
			nmesh->insert_vtx(_msh->GetVerticesEditable()[oriid].pos);
			idcounter++;
			idx_mapping.insert(std::make_pair(std::make_pair(oriid, i), idcounter));
		}
	}

	std::vector<std::vector<int>> cc_contained_tris;
	cc_contained_tris.resize(_gh.node_list.size(), std::vector<int>{});
	for (int i = 0; i < _msh->GetNumTriangles(); i++) {
		int ccid = _msh->GetTriangles()[i].component_id;
		cc_contained_tris[ccid].push_back(i);
	}

	auto get_new_id = [](int vtx_id, int faceid, std::map<int, std::vector<std::vector<int>>> vertex_split_groups, std::map<std::pair<int, int>, int> &idx_mapping) {
		auto it = vertex_split_groups.find(vtx_id);
		if (it == vertex_split_groups.end())
			return vtx_id;

		int groupid = -1;
		for (int i = 0; i < it->second.size(); i++) {
			for (int j = 0; j < it->second[i].size(); j++) {
				if (it->second[i][j] == faceid) {
					groupid = i;
					break;
				}
			}
			if (groupid != -1) break;
		}

		return idx_mapping[std::make_pair(vtx_id, groupid)];
	};

	for (int j = 0; j < ncclist.size(); j++) {

		for (int m = 0; m < ncclist[j].size(); m++) {
			int ccid = ncclist[j][m];
#ifdef _OUTPUT_DEBUG_LOGS
			if (ccdir[ccid] == -1) {
				std::cout << "Error";
#if defined(_RECORD_ERROR_LOGS)
				(*_log_off) << _file_under_processing << " ERROR write_mesh_change_topo_with_solution\n";
				(*_log_off).flush();
#endif
			}
#endif

			bool positive = ccdir[ccid] == 0;

			for (int m = 0; m < cc_contained_tris[ccid].size(); m++) {
				int triid = cc_contained_tris[ccid][m];
				int ori_ccid = _msh->GetTriangles()[triid].component_id;
				if (del_cc_ls[ori_ccid]) continue;
				if (positive) {
					nmesh->insert_face(std::vector<int>{get_new_id(_msh->GetTriangles()[triid].vertices[0], triid, vertex_split_groups, idx_mapping),
						get_new_id(_msh->GetTriangles()[triid].vertices[1], triid, vertex_split_groups, idx_mapping),
						get_new_id(_msh->GetTriangles()[triid].vertices[2], triid, vertex_split_groups, idx_mapping)});
				}
				else {
					nmesh->insert_face(std::vector<int>{get_new_id(_msh->GetTriangles()[triid].vertices[2], triid, vertex_split_groups, idx_mapping),
						get_new_id(_msh->GetTriangles()[triid].vertices[1], triid, vertex_split_groups, idx_mapping),
						get_new_id(_msh->GetTriangles()[triid].vertices[0], triid, vertex_split_groups, idx_mapping)});
				}
			}
		}
	}

	delete _msh;
	_msh = nmesh;
	_msh->update_mesh_properties();
	_msh->mark_component_with_coherence();
	return;
}

bool Optimizer::postProcess(std::vector<int>& x, POSTPRO_TYPE post_type,
	std::vector<bool> & del_cc_ls, std::vector<std::vector<std::vector<int>>> &ccid_footprint)
{

	//true- no cc deleted
	std::vector<int> ccdir;
	std::vector<std::vector<int>> ncclist;
	get_directioend_component_list(x, ccdir, ncclist);

	ccid_footprint.push_back(ncclist);

	std::vector<int> oldcc_to_newcc;
	oldcc_to_newcc.resize(_gh.node_list.size(), -1);

	//build map from old cc to new cc
	for (int i = 0; i < ncclist.size(); i++) {
		for (int j = 0; j < ncclist[i].size(); j++) {
			oldcc_to_newcc[ncclist[i][j]] = i;
		}
	}

	VisualProcesser vp;
	vp.setMesh(_msh);

	std::vector<int> posPx, negPx;
	vp.get_unary_score_all_campos(posPx, negPx);

	std::vector<int> ccposPx, ccnegPx;
	ccposPx.resize(_gh.node_list.size(), 0);
	ccnegPx.resize(_gh.node_list.size(), 0);

	for (int i = 0; i < posPx.size(); i++) {
		int ccid = _msh->GetTriangles()[i].component_id;
		ccposPx[ccid] += posPx[i];
		ccnegPx[ccid] += negPx[i];
	}

	//std::vector<int> del_cc_ls;
	del_cc_ls.clear();
	del_cc_ls.resize(ccposPx.size(), false);
	bool reconstruct_flag = false;
	if (post_type == DELETE_ALL_IN) {
		for (int i = 0; i < ccposPx.size(); i++) {
			if (ccposPx[i] <= _TREAT_PIX_SMALLER_THAN_AS_ZERO && ccnegPx[i] <= _TREAT_PIX_SMALLER_THAN_AS_ZERO) {
				del_cc_ls[i] = true;
				reconstruct_flag = true;
			}
		}
	}
	else {
		for (int i = 0; i < ccposPx.size(); i++) {
			if (ccposPx[i] <= _TREAT_PIX_SMALLER_THAN_AS_ZERO && ccnegPx[i] <= _TREAT_PIX_SMALLER_THAN_AS_ZERO) {
				bool push_flag = false;
				for (int j = 0; j < ncclist[oldcc_to_newcc[i]].size(); j++) {
					int thiscc = ncclist[oldcc_to_newcc[i]][j];
					if (thiscc != i &&
						((ccposPx[thiscc] <= _TREAT_PIX_SMALLER_THAN_AS_ZERO && ccnegPx[thiscc] > _TREAT_PIX_SMALLER_THAN_AS_ZERO)
							|| (ccposPx[thiscc] > _TREAT_PIX_SMALLER_THAN_AS_ZERO && ccnegPx[thiscc] <= _TREAT_PIX_SMALLER_THAN_AS_ZERO))) {
						push_flag = true;
						break;
					}
				}
				if (push_flag) {
					del_cc_ls[i] = true;
					reconstruct_flag = true;
				}
			}
		}
	}

	//delete faces but donot change topo
	if (reconstruct_flag) {
		BlackMesh::BlackMesh<double> *nmesh = new BlackMesh::BlackMesh<double>();


		nmesh->copy_vertice_from(_msh);

#ifdef _OUTPUT_OPT_DEL_AS_MESH
		BlackMesh::BlackMesh delepart;
		delepart.copy_vertice_from(_msh);
#endif

		for (int i = 0; i < _msh->GetNumTriangles(); i++) {
			int ccid = _msh->GetTriangles()[i].component_id;
			if (!del_cc_ls[ccid]) {
				nmesh->insert_face(_msh->GetTrianglesEditable()[i].vertices);
			}
#ifdef _OUTPUT_OPT_DEL_AS_MESH
			else {
				delepart.insert_face(_msh->GetTrianglesEditable()[i].vertices);
			}
#endif
		}

#ifdef _OUTPUT_OPT_DEL_AS_MESH
		MeshIO::writeOBJ("_debug_opt_deleted.obj", delepart);
#endif

		delete _msh;
		_msh = nmesh;
		_msh->update_mesh_properties();
		_msh->mark_component_with_coherence();

		//MeshIO::writeOBJ("ttt.obj", (*_msh));

		return !reconstruct_flag;
	}
	else {
		//MeshIO::writeOBJ("ttt.obj", (*_msh));
		return !reconstruct_flag;
	}
}



void Optimizer::write_mesh_grouped_with_solution(std::vector<int> &x) {
	std::vector<int> ccdir;
	std::vector<vector<int>> ncclist;
	get_directioend_component_list(x, ccdir, ncclist);

	ofstream off("_debug_output.obj");
	off << "mtllib output.mtl\n";
	ofstream offm("output.mtl");
	for (int i = 0; i < _msh->GetNumVertices(); i++) {
		off << "v " << _msh->GetVertices()[i].pos[0]
			<< " " << _msh->GetVertices()[i].pos[1]
			<< " " << _msh->GetVertices()[i].pos[2] << "\n";
	}

	std::vector<std::vector<int>> cc_contained_tris;
	cc_contained_tris.resize(_gh.node_list.size(), std::vector<int>{});
	for (int i = 0; i < _msh->GetNumTriangles(); i++) {
		int ccid = _msh->GetTriangles()[i].component_id;
		cc_contained_tris[ccid].push_back(i);
	}

	for (int j = 0; j < ncclist.size(); j++) {
		off << "g cc" << j << "\n";
		off << "usemtl mcc" << j << "\n";

		for (int m = 0; m < ncclist[j].size(); m++) {
			int ccid = ncclist[j][m];
#ifdef _OUTPUT_DEBUG_LOGS
			if (ccdir[ccid] == -1) {
				std::cout << "Error";
#if defined(_RECORD_ERROR_LOGS)
				(*_log_off) << _file_under_processing << " ERROR write_mesh_grouped_with_solution\n";
				(*_log_off).flush();
#endif
#endif
			}
			bool positive = ccdir[ccid] == 0;

			for (int m = 0; m < cc_contained_tris[ccid].size(); m++) {
				int triid = cc_contained_tris[ccid][m];
				if (positive) {
					off << "f " << _msh->GetTriangles()[triid].vertices[0] + 1 << " "
						<< _msh->GetTriangles()[triid].vertices[1] + 1 << " "
						<< _msh->GetTriangles()[triid].vertices[2] + 1 << " \n";
				}
				else {
					off << "f " << _msh->GetTriangles()[triid].vertices[2] + 1 << " "
						<< _msh->GetTriangles()[triid].vertices[1] + 1 << " "
						<< _msh->GetTriangles()[triid].vertices[0] + 1 << " \n";
				}
			}
		}

		offm << "newmtl mcc" << j << "\n";
		offm << "\tNs 10.0\n\tNi 1.5\n\td 1.0\n\tTr 0.0\n\tillum 2\n\tKa 1.0 1.0 1.0\n\tKd "
			<< color_table[j % 17][0] << " " << color_table[j % 17][1] << " " << color_table[j % 17][2] << "\n\tKs 0.0 0.0 0.0\n\tKe 0.0 0.0 0.0\n";
	}
	off.close();
	offm.close();
}

void Optimizer::write_mesh_change_topo_with_solution(std::vector<int>& x)
{
	//x is edge-based
	auto build_oriedge_to_graphedge_map = [](map<int, int> &res, BlackMesh::BlackMesh<double> *msh, Primal_Dual_graph *gh) {
		res.clear();

		for (int i = 0; i < gh->edge_list.size(); i++) {

			res.insert(std::make_pair(gh->edge_list[i].original_edge_id, i));

		}

#ifdef _OUTPUT_DEBUG_LOGS
		for (int i = 0; i < gh->edge_list.size(); i++) {
			if (res[gh->edge_list[i].original_edge_id] != i) {
				std::cout << "The map between edges is not matched\n";
#if defined(_RECORD_ERROR_LOGS)
				(*_log_off) << _file_under_processing << " ERROR write_mesh_change_topo_with_solution\n";
				(*_log_off).flush();
#endif
			}
		}
#endif
	};

	auto get_incident_non_manifold_edges = [](map<int, int> &mp, std::vector<int> &ori_edgeid, std::vector<int> &gh_edgeid) {
		gh_edgeid.clear();
		for (int i = 0; i < ori_edgeid.size(); i++) {
			auto it = mp.find(ori_edgeid[i]);
			if (it == mp.end()) {
				//do nothing
			}
			else {
				gh_edgeid.push_back(it->second);
			}
		}
	};

	auto get_eq_face_set_on_vtx = [&](int vertexid, const map<int, int> &msh2ghmap, const std::vector<int>& x, std::vector<std::vector<int>> &output_groups) {
		std::set<int> visited_face;
		std::set<int> all_face;
		//get all faces incident to vtx
		auto incident_edges_to_vtx = _msh->GetVertices()[vertexid].incident_edges;
		for (int i = 0; i < incident_edges_to_vtx.size(); i++) {
			all_face.insert(_msh->GetEdges()[incident_edges_to_vtx[i]].incident_faces.begin(),
				_msh->GetEdges()[incident_edges_to_vtx[i]].incident_faces.end());
		}


		//go through all faces
		output_groups.clear();
		for (auto it = all_face.begin(); it != all_face.end(); it++) {
			int faceid = *it;
			if (visited_face.find(faceid) != visited_face.end())continue;

			std::queue<int> list;
			list.push(faceid);
			visited_face.insert(faceid);
			output_groups.push_back(std::vector<int>{faceid});
			int siz = output_groups.size();
			int groupid = siz - 1;

			while (list.size() != 0) {
				int current_face = list.front();
				list.pop();

				for (int i = 0; i < _msh->GetTriangles()[current_face].incident_edges.size(); i++) {
					int edgeidmsh = _msh->GetTriangles()[current_face].incident_edges[i];

					if (_msh->GetEdges()[edgeidmsh].is_nonmanifold_edge == false) {
						//not nedge
						for (int j = 0; j < _msh->GetEdges()[edgeidmsh].incident_faces.size(); j++) {
							int query_face_id = _msh->GetEdges()[edgeidmsh].incident_faces[j];
							if (visited_face.find(query_face_id) != visited_face.end())continue;
							//if query face have vtx
							bool hit = false;

							if (all_face.find(query_face_id) != all_face.end())hit = true;

							if (hit) {
								list.push(query_face_id);
								visited_face.insert(query_face_id);
								output_groups[groupid].push_back(query_face_id);
							}
						}

					}
					else {
						//nedge
						int edgeidgh = msh2ghmap.find(edgeidmsh)->second;
						int degree = _gh.edge_list[edgeidgh].incident_components.size();
						auto pickedconf = configuration_list_all_degrees[degree - 2][x[edgeidgh]];

						//get the incident id of current face
						int incident_id_of_current_face = -1;
						for (int j = 0; j < _msh->GetEdges()[edgeidmsh].incident_faces.size(); j++) {
							if (_msh->GetEdges()[edgeidmsh].incident_faces[j] == current_face) {
								incident_id_of_current_face = j;
								break;
							}
						}

						//find the dual edge picked that have such incident id
						int anotherId = -2;
						for (int j = 0; j < pickedconf.size(); j++) {
							auto picked_dualedge = dual_edge_list_all_degrees[degree - 2][pickedconf[j]];

							if (picked_dualedge.fi == incident_id_of_current_face) {
								if (anotherId != -2) {
									std::cout << "ERROR HERR\n";
#if defined(_RECORD_ERROR_LOGS)
									(*_log_off) << _file_under_processing << " ERROR write_mesh_change_topo_with_solution\n";
									(*_log_off).flush();
#endif
								}
								anotherId = picked_dualedge.fj;
							}
							else if (picked_dualedge.fj == incident_id_of_current_face) {
								if (anotherId != -2) {
									std::cout << "ERROR HERR\n";
#if defined(_RECORD_ERROR_LOGS) 
									(*_log_off) << _file_under_processing << " ERROR write_mesh_change_topo_with_solution\n";
									(*_log_off).flush();
#endif
								}
								anotherId = picked_dualedge.fi;
							}
						}

						if (anotherId == -2) {
							std::cout << "ERROR HERR\n";
#if defined(_RECORD_ERROR_LOGS)
							(*_log_off) << _file_under_processing << " ERROR write_mesh_change_topo_with_solution\n";
							(*_log_off).flush();
#endif
						}
						else if (anotherId == -1) {
							//current face is break to the nedge
						}
						else {
							int query_face_id = _msh->GetEdges()[edgeidmsh].incident_faces[anotherId];
							if (visited_face.find(query_face_id) != visited_face.end())continue;
							//if query face have vtx
							bool hit = false;
							if (all_face.find(query_face_id) != all_face.end())hit = true;
							if (hit) {
								list.push(query_face_id);
								visited_face.insert(query_face_id);
								output_groups[groupid].push_back(query_face_id);
							}
						}

					}

				}

			}
		}

	};


	std::map<int, int> old_to_new_edge;
	build_oriedge_to_graphedge_map(old_to_new_edge, _msh, &_gh);

	std::vector<std::vector<std::vector<int>>> vertex_split_groups; //vid as index
	vertex_split_groups.resize( _msh->GetNumVertices(), std::vector<std::vector<int>>{});
	for (int i = 0; i < _msh->GetNumVertices(); i++) {
		get_eq_face_set_on_vtx(i, old_to_new_edge, x, vertex_split_groups[i]);
	}



	/*****/
		//get component list
	std::vector<int> ccdir;
	std::vector<std::vector<int>> ncclist;
	get_directioend_component_list(x, ccdir, ncclist);

	const char *objsuf = ".obj";
	const char *mtlsuf = ".mtl";
	const size_t lenobj = Folder_path.length() + Output_filename.length() + strlen(objsuf);
	const size_t lenmtl = Folder_path.length() + Output_filename.length() + strlen(mtlsuf);

	char *n_strobj = new char[lenobj + 1];
	strcpy(n_strobj, Folder_path.c_str());
	strcat(n_strobj, Output_filename.c_str());
	strcat(n_strobj, objsuf);

	string mtlname = "commonMTL";
	char *n_strmtl = new char[lenmtl + 1];
	strcpy(n_strmtl, Folder_path.c_str());
	strcat(n_strmtl, mtlname.c_str());
	strcat(n_strmtl, mtlsuf);

	ofstream off(n_strobj);
	off << "mtllib  commonMTL.mtl\n";
	ofstream offm(n_strmtl);
	for (int i = 0; i < _msh->GetNumVertices(); i++) {
		off << std::setprecision(18) << "v " << _msh->GetVertices()[i].pos[0]
			<< " " << _msh->GetVertices()[i].pos[1]
			<< " " << _msh->GetVertices()[i].pos[2] << "\n";
	}

	//write extra vtx
	std::map<std::pair<int, int>, int> idx_mapping;//ori_vtxid,ccid-newvtxid
	int idcounter = _msh->GetNumVertices();
	idcounter--;

	for (int ii = 0; ii < _msh->GetNumVertices(); ii++) {
		//0th, push original id


		int oriid = ii;
		int groupid = 0;
		idx_mapping.insert(std::make_pair(std::make_pair(oriid, 0), oriid));


		for (int i = 1; i < vertex_split_groups[ii].size(); i++) {
			int oriid = ii;
			off << std::setprecision(18) << "v " << _msh->GetVertices()[oriid].pos[0]
				<< " " << _msh->GetVertices()[oriid].pos[1]
				<< " " << _msh->GetVertices()[oriid].pos[2] << "\n";
			idcounter++;


			idx_mapping.insert(std::make_pair(std::make_pair(oriid, i), idcounter));


		}
	}

	std::vector<std::vector<int>> cc_contained_tris;
	cc_contained_tris.resize(_gh.node_list.size(), std::vector<int>{});
	for (int i = 0; i < _msh->GetNumTriangles(); i++) {
		int ccid = _msh->GetTriangles()[i].component_id;
		cc_contained_tris[ccid].push_back(i);
	}

	struct what_is_nid {
		bool need_split;
		std::map<int, int> fid_to_nid;
		what_is_nid() :need_split(false) {

		};
	};

	auto get_new_id = [](int vid, int fid, std::vector<what_is_nid> &nidlist) {
		if (nidlist[vid].need_split) {
			auto it = nidlist[vid].fid_to_nid.find(fid);
			if (it == nidlist[vid].fid_to_nid.end()) {
				assert(false);
				return -1;
			}
			else {
				return it->second;
			}
		}
		else {
			return vid;
		}
	};

	//construct o(1) visit for nid


	std::vector<what_is_nid> nidlist;
	nidlist.resize(_msh->GetNumVertices());

	for (int ii = 0; ii < _msh->GetNumVertices(); ii++) {
		int vid = ii;
		nidlist[vid].need_split = true;

		for (int i = 0; i < vertex_split_groups[ii].size(); i++) {
			int nid = idx_mapping[std::make_pair(vid, i)];

			for (int j = 0; j <vertex_split_groups[ii][i].size(); j++) {
				int fid = vertex_split_groups[ii][i][j];
				nidlist[vid].fid_to_nid.insert(std::make_pair(fid, nid));
			}

		}
	}

	for (int j = 0; j < ncclist.size(); j++) {
		off << "g cc" << j << "\n";
		off << "usemtl mcc" << j % 17 << "\n";

		for (int m = 0; m < ncclist[j].size(); m++) {
			int ccid = ncclist[j][m];
#ifdef _OUTPUT_DEBUG_LOGS
			if (ccdir[ccid] == -1) {
				std::cout << "Error";
#if defined(_RECORD_ERROR_LOGS)
				(*_log_off) << _file_under_processing << " ERROR write_mesh_change_topo_with_solution\n";
				(*_log_off).flush();
#endif
			}
#endif

			bool positive = ccdir[ccid] == 0;

			for (int m = 0; m < cc_contained_tris[ccid].size(); m++) {
				int triid = cc_contained_tris[ccid][m];
				if (positive) {
					off << "f " << get_new_id(_msh->GetTriangles()[triid].vertices[0], triid, nidlist) + 1 << " "
						<< get_new_id(_msh->GetTriangles()[triid].vertices[1], triid, nidlist) + 1 << " "
						<< get_new_id(_msh->GetTriangles()[triid].vertices[2], triid, nidlist) + 1 << " \n";
				}
				else {
					off << "f " << get_new_id(_msh->GetTriangles()[triid].vertices[2], triid, nidlist) + 1 << " "
						<< get_new_id(_msh->GetTriangles()[triid].vertices[1], triid, nidlist) + 1 << " "
						<< get_new_id(_msh->GetTriangles()[triid].vertices[0], triid, nidlist) + 1 << " \n";
				}
			}
		}


	}

	for (int j = 0; j < 17; j++) {
			offm << "newmtl mcc" << j % 17 << "\n";
			offm << "\tNs 10.0\n\tNi 1.5\n\td 1.0\n\tTr 0.0\n\tillum 2\n\tKa 1.0 1.0 1.0\n\tKd "
				<< color_table[j % 17][0] << " " << color_table[j % 17][1] << " " << color_table[j % 17][2] << "\n\tKs 0.900 0.900 0.900\n\tKe 0.0 0.0 0.0\n";

	}

	off.close();
	offm.close();
	/**********/

	delete[]n_strmtl;
	delete[]n_strobj;
	return;
}


bool Optimizer::check_global_laplacian()
{
	return false;
}

bool Optimizer::check_ori_cons_with_confid()
{
	//std::vector<map<int, map<int, std::vector<orientation_Expr_with_id>>>> ori_cons_with_confid;
	//ccid,edge id, incident id, cons// the incident id is to seperate the case for degenerate Lap
	//the direction is based on the original normal of cc. dir true mean keep its original normal
	//the direction of dual edge is based on the smaller incident id. dir=0 means keep the smaller one. dir =1 means flip the smaller one;
	//connected face normal will be coherent following manifold property

	for (int ccid = 0; ccid < ori_cons_with_confid.size(); ccid++) {
		for (auto edgeit = ori_cons_with_confid[ccid].begin(); edgeit != ori_cons_with_confid[ccid].end(); edgeit++) {
			int curveid = edgeit->first;
			for (auto infaceit = edgeit->second.begin(); infaceit != edgeit->second.end(); infaceit++) {
				int incidentid = infaceit->first;
				for (auto confit = infaceit->second.begin(); confit != infaceit->second.end(); confit++) {
					int pushed_curveid = confit->edgeid;
					int pushed_conf = confit->confid;
					bool induced_face_dir_is_pos = confit->positive_flg;

					if (curveid != pushed_curveid) {
						std::cout << "Stored edge id is not consistent\n";
#if defined(_RECORD_ERROR_LOGS) 
						(*_log_off) << _file_under_processing << " ERROR Stored edge id is not consistent\n";
						(*_log_off).flush();
#endif
					}


					//check the conf of the edge surely connect the cc
					auto eh = _gh.edge_list[_gh.non_manifold_curves[curveid][0]];
					int degree = eh.incident_components.size();
					auto picked_dual_edges = configuration_list_all_degrees[degree - 2][pushed_conf];


					bool conf_do_connect_ccid = false;
					bool cc_dir_correct = false;

					for (int id = 0; id < picked_dual_edges.size(); id++) {
						int picked_dual = picked_dual_edges[id];
						int fi = dual_edge_list_all_degrees[degree - 2][picked_dual].fi;
						int fj = dual_edge_list_all_degrees[degree - 2][picked_dual].fj;
						int dir = dual_edge_list_all_degrees[degree - 2][picked_dual].dir;

						//the incident id fi should less than fj by construction
						if (fj != -1) {
							assert(fi < fj);

							int fidi = eh.incident_components[fi];
							int fidj = eh.incident_components[fj];

							if ((fidi == ccid&&fi == incidentid) || (fidj == ccid&&fj == incidentid)) {
								conf_do_connect_ccid = true;
							}
							else {
								continue;
							}

							if (ccid == fidi&&fi == incidentid) {
								//the edge dir is consistent with face
								if ((dir == 0 && induced_face_dir_is_pos) || (dir == 1 && !induced_face_dir_is_pos)) {
									//pass
									cc_dir_correct = true;
								}
								else {
									
								}
							}

							if (ccid == fidj&&fj == incidentid) {
								//the direction of cc should rely on the dir of fidi
								int orif_fi = _msh->GetEdges()[eh.original_edge_id].incident_faces[fi];
								int orif_fj = _msh->GetEdges()[eh.original_edge_id].incident_faces[fj];

								int flg = if_face_coherant_based_on_edge(_msh->GetEdges()[eh.original_edge_id].vertices.first, _msh->GetEdges()[eh.original_edge_id].vertices.second,
									_msh->GetTriangles()[orif_fi], _msh->GetTriangles()[orif_fj]);

								if (flg == 1) {
									if ((dir == 0 && induced_face_dir_is_pos) || (dir == 1 && !induced_face_dir_is_pos)) {
										//pass
										cc_dir_correct = true;
									}
									else {
										
									}
								}

								if (flg == 0) {
									if ((dir == 1 && induced_face_dir_is_pos) || (dir == 0 && !induced_face_dir_is_pos)) {
										//pass
										cc_dir_correct = true;

									}
									else {
										
									}
								}
							}

						}
						else {
							int fidi = eh.incident_components[fi];
							if (fidi == ccid&&fi == incidentid) {
								conf_do_connect_ccid = true;
							}
							else {
								continue;
							}
							//
							if ((dir == 0 && induced_face_dir_is_pos) || (dir == 1 && !induced_face_dir_is_pos)) {
								//pass
								cc_dir_correct = true;
							}
							else {
								
							}
						}
					}

					if (conf_do_connect_ccid == false) {
						cout << "The conf do not connect cc\n";
#if defined(_RECORD_ERROR_LOGS)
						(*_log_off) << _file_under_processing << " ERROR The cc dir is not correct\n";
						(*_log_off).flush();
#endif
					}
					if (cc_dir_correct == false) {
						cout << "The cc dir is not correct\n";
#if defined(_RECORD_ERROR_LOGS)
						(*_log_off) << _file_under_processing << " ERROR The cc dir is not correct\n";
						(*_log_off).flush();
#endif
					}

				}
			}
		}
	}

	return false;

}

void Optimizer::show_solution(std::vector<std::vector<double>>& x)
{

	std::cout << "*****************Solution Check*********************\n";
	std::cout << "CurvId\t\tEdge0\t\tIncidentCC\t\t\t\tConf\t\tConfDetail\n";

	for (int i = 0; i < x.size(); i++) {
		int curveid = i;
		int edgeid = _gh.non_manifold_curves[i][0];
		auto incidentCC = _gh.edge_list[edgeid].incident_components;
		int Conf = -1;
		for (int j = 0; j < x[i].size(); j++) {
			if (abs(x[i][j] - 1.0) < 1e-7) {
				Conf = j;
				break;
			}
		}
		int degree = incidentCC.size();
		auto confDualEdgesId = configuration_list_all_degrees[degree - 2][Conf];

		std::cout << curveid << "\t\t";
		std::cout << edgeid << "\t\t";
		for (int j = 0; j < incidentCC.size(); j++) {
			std::cout << incidentCC[j] << " ";
		}
		std::cout << "\t\t\t\t";

		std::cout << Conf << "\t\t";

		for (int j = 0; j < confDualEdgesId.size(); j++) {
			int dualid = confDualEdgesId[j];
			std::cout << "(" << dual_edge_list_all_degrees[degree - 2][dualid].fi << "," <<
				dual_edge_list_all_degrees[degree - 2][dualid].fj << "," <<
				dual_edge_list_all_degrees[degree - 2][dualid].dir << ")";
		}
		std::cout << "\n";
	}


}



