#pragma once
#include "BlackMesh.h"
#include "Eigen/Dense"

#include <cmath>
#include <cassert>
#include <functional>
#include "SparseTensor.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <set>
#include "Mosek_c_wrapper.h"
#include "Primal_Dual_graph.h"
#include "MeshIO.h"
#include "MRF_solver.h"
#include "Tribunal.h"

class Optimizer
{
private:
	struct dual_edge_st {
		int fi;
		int fj;
		int dir;//0-positive
		dual_edge_st(int fi_, int fj_, int dir_) :fi(fi_), fj(fj_), dir(dir_) {};
	};

	struct orientation_Expr_with_id {
		bool positive_flg;//0- negative
		int edgeid;
		int confid;
		orientation_Expr_with_id(bool positive_flg_, int edgeid_, int confid_) :positive_flg(positive_flg_), edgeid(edgeid_), confid(confid_) {}
	};

	BlackMesh::BlackMesh<double> *_msh;
	Primal_Dual_graph _gh;
	 char* Output_filename;
	 char* Folder_path;

	 int TIME_MIN;
	 int TIME_MAX;


	std::vector<std::vector<std::vector<int>>> configuration_list_all_degrees;//degree-configuration id-store which dual edge is picked
	std::vector<std::vector<dual_edge_st>> dual_edge_list_all_degrees;//degree-dual edge id

	std::vector<map<int, map<int,std::vector<orientation_Expr_with_id>>>> ori_cons_with_confid;//ccid,edge id, incident id, cons// the incident id is to seperate the case for degenerate Lap

	std::map<std::pair<int, int>, std::vector<Tri>> global_laplacian;//ccidi, ccidj-> edgeid, confid, val//ccidi<=ccidj triangle matrix

	Mosek_c_wrapper base_msk_solver;
	map<std::pair<int, int>, int> var_map;

	int n_opt_vars;

	std::vector<SpMat> trans_conf2dual_eign;

	
	int Optimizer::if_face_coherant_based_on_edge(int vid_start, int vid_end,
		const BlackMesh::BlackMesh<double>::Triangle &tr1, const BlackMesh::BlackMesh<double>::Triangle &tr2);

	//init
	void generate_transformation_matrices();
	int  Optimizer::get_dual_edge_id(int fi, int fj, int degree, int direction);
	void  Optimizer::add_one_dual_edge_on_current_configurations(std::vector<std::vector<int>> &confs_res, std::vector<dual_edge_st> &dual_edge_list);


	//for score assemble and rounding
	//std::vector<std::pair<double, std::vector<std::pair<int,int>>>> unary_score;// 2*ccid+0/1 -> score val, vec(edgeid, confid)
	//std::map<int, vector<Binary_element>> unary_score_from_binary;//2*ccid+0/1 -> vec(edgeid, confid, val)
	//std::map<std::pair<int, int>, vector<vector<Binary_element>>> binary_score;// ccidi, ccidj -> vec(00/01/10/11)(edgeid, confid, val)
	//std::map<std::pair<int, int>, vector<std::pair<int,int>>> dual_induced_dir;// ccidi, ccidj

	std::map<std::pair<int, int>, Binary_Element> binary_score;
	std::map<std::pair<int, int>, Binary_Element> unary_score_from_binary;
	std::vector<Unary_Element> unary_score;

	inline void insert_to_binary_score(int curveid, int confid, int dualid, int inidi, int inidj, int diri, int dirj) {
		auto it = binary_score.insert(std::make_pair(std::make_pair(curveid, dualid), Binary_Element()));
		if (it.second == true) {
			it.first->second.inidi = inidi;
			it.first->second.inidj = inidj;
			it.first->second.diri = diri;
			it.first->second.dirj = dirj;
			it.first->second.confid.insert(confid);
		}
		else {
			it.first->second.confid.insert(confid);
#ifdef _OUTPUT_DEBUG_LOGS
			if (it.first->second.inidi != inidi) {
				std::cout << "ERROR\n";
#if defined(_RECORD_ERROR_LOGS)
				(*_log_off) << _file_under_processing << " ERROR 1\n";
				(*_log_off).flush();
#endif
			}
			if (it.first->second.inidj != inidj) {
				std::cout << "ERROR\n";
#if defined(_RECORD_ERROR_LOGS)
				(*_log_off) << _file_under_processing << " ERROR 1\n";
				(*_log_off).flush();
#endif
			}
			if (it.first->second.diri != diri) {
				std::cout << "ERROR\n";
#if defined(_RECORD_ERROR_LOGS) 
				(*_log_off) << _file_under_processing << " ERROR 1\n";
				(*_log_off).flush();
#endif
			}
			if (it.first->second.dirj != dirj) {
				std::cout << "ERROR\n";
#if defined(_RECORD_ERROR_LOGS)
				(*_log_off) << _file_under_processing << " ERROR 1\n";
				(*_log_off).flush();
#endif
			}
#endif
		}
	}
	inline void insert_to_unary_score(int curvid, int confid, int ccid, int dir) {
		if (dir == 0) {
			unary_score[ccid].confid_pos.insert(std::make_pair(curvid, confid));
		}
		else {
			unary_score[ccid].confid_neg.insert(std::make_pair(curvid, confid));
		}
	}
	inline void insert_to_unary_from_binary_score(int curveid, int confid, int dualid, int inidi, int inidj, int dir) {
		auto it = unary_score_from_binary.insert(std::make_pair(std::make_pair(curveid, dualid), Binary_Element()));
		if (it.second == true) {
			it.first->second.inidi = inidi;
			it.first->second.inidj = inidj;
			it.first->second.diri = dir;
			it.first->second.dirj = -1;
			it.first->second.confid.insert(confid);
		}
		else {
			it.first->second.confid.insert(confid);
#ifdef _OUTPUT_DEBUG_LOGS
			if (it.first->second.inidi != inidi) {
				std::cout << "ERROR\n";
#if defined(_RECORD_ERROR_LOGS)
				(*_log_off) << _file_under_processing << " ERROR 1\n";
				(*_log_off).flush();
#endif
			}
			if (it.first->second.inidj != inidj) {
				std::cout << "ERROR\n";
#if defined(_RECORD_ERROR_LOGS) 
				(*_log_off) << _file_under_processing << " ERROR 1\n";
				(*_log_off).flush();
#endif
			}
			if (it.first->second.diri != dir) {
				std::cout << "ERROR\n";
#if defined(_RECORD_ERROR_LOGS)
				(*_log_off) << _file_under_processing << " ERROR 1\n";
				(*_log_off).flush();
#endif
			}
			if (it.first->second.dirj != -1) {
				std::cout << "ERROR\n";
#if defined(_RECORD_ERROR_LOGS) 
				(*_log_off) << _file_under_processing << " ERROR 1\n";
				(*_log_off).flush();
#endif
			}
#endif
		}
	}

	void solve_ccdir_on_current_sol(std::vector<vector<double>> &sol, std::map<int, int> *node, std::vector<int> &ccdir);
	int get_max_feasiable_confid(std::vector<vector<double>> &sol, std::vector<int> &ccdir, int edgeid);
	void MRF_rounding(std::vector<vector<double>> &sol, std::map<int, int> *node, std::vector<vector<double>> &rounded);

	std::vector<std::pair<int,int>> obj_vid;
	std::vector<double> obj_val;

	std::map<std::pair<int, int>, double> valid_val_unary;
	std::map<std::pair<int, int>, double> valid_val_binary;

	//
	MeshSampling _sampler;
	std::vector<double> non_manifold_curve_length;
	inline void update_sampler(int ave_sampling = 1000, int min_sampling = 20) {
		_sampler.setMesh(_msh);
		_sampler.generate_sampling_points(ave_sampling, min_sampling);
	};
	void update_curve_length();

	double _sWeight;
	double _uWeight;
	double _bWeight;

	double total_lap_weight;

public:

	//init
	inline Primal_Dual_graph* getPrimalDualGraph() { return &_gh; }
	inline void setMesh(BlackMesh::BlackMesh<double> *msh) { _msh = msh; }
	void init();
	void split_unorientated_cases(std::vector<std::vector<int>> &confs_res, std::vector<dual_edge_st> &dual_edge_list, int degree, int dual_siz, std::vector<std::vector<int>>& output);

	//stack edges to curves
	void stack_non_manifold_curves(std::vector<std::vector<int>> &nonmanifold_curves);


	//others
	void turn_conf_to_stacked(std::vector<std::vector<double>> &xin, std::vector<int> &xout);////

	inline void setOutputName(  char* folder,  char* name) { Folder_path = folder; Output_filename = name; }

	inline void setTimer(int min, int max) {
		TIME_MIN = min;
		TIME_MAX = max;
	}

	//assemnble and solve
	inline int pre_assemble() {

		_gh.build_connection_graph(_msh);
		if (_gh.edge_list.size() == 0)return 1;
		//_gh.update_combined_list();
		init();
		assemble_contraints();

		build_variable_idx_map(var_map);
		assemble_solver(base_msk_solver, var_map);

		return 0;
	}
	inline int pre_assemble_curve_based() {


		_gh.clear();
		_msh->clearNumComponent();
		_gh.build_connection_graph(_msh);
		if (_gh.edge_list.size() == 0) {
				
			return 1;
		}


		stack_non_manifold_curves(_gh.non_manifold_curves);
		std::vector<int> edgeidNeedToSplit;
		for (int i = 0; i < _gh.non_manifold_curves.size(); i++) {
			if (_gh.non_manifold_curves[i].size() == 1) {
				edgeidNeedToSplit.push_back(_gh.edge_list[_gh.non_manifold_curves[i][0]].original_edge_id);
			}
		}

		
		if (edgeidNeedToSplit.size() > 0) {
			_msh->split_edge(edgeidNeedToSplit);
			_gh.clear();
			_gh.build_connection_graph(_msh);
			stack_non_manifold_curves(_gh.non_manifold_curves);
			std::cout << edgeidNeedToSplit.size()<<" Edges Splited\n";
		}



		init();
		update_sampler();

		update_curve_length();

		assemble_contraints();
		update_unary_and_binary_score();
		build_variable_idx_map_curve_based(var_map);
		std::cout << "There are " << var_map.size() << " variables\n";

		assemble_imcompatible_confs();

		assemble_solver(base_msk_solver, var_map);

		return 0;
	}
	void Optimizer::assemble_orientation_cons(int edgeid);
	void Optimizer::assemble_orientation_cons_curve_based(int curveid);
	void Optimizer::assemble_contraints_on_one_nonmanifold_edge(int i);
	void Optimizer::assemble_contraints();

	std::vector<std::pair<int, int>> imcompatible_confs;//curvid confid
	void Optimizer::assemble_imcompatible_confs();
	void Optimizer::append_imcp_confs(map<std::pair<int, int>, int>& mp, Mosek_c_wrapper &solver,
		std::vector<std::pair<int, int>> &imcompatible_confs);

	void append_cons_on_all_edges(std::map<int, int>& node);
	
	inline void assemble_solver(Mosek_c_wrapper &solver, map<std::pair<int, int>, int> &idx_mp) {
		solver.make_task();

		//build idx mapping

		if (_gh.node_list.size() > 1) {

			solver.append_vars(n_opt_vars, 0);

			append_orientation_cons_curve_based(idx_mp, solver);

			append_x_regular_cons_curve_based(idx_mp, solver);

			append_imcp_confs(idx_mp, solver, imcompatible_confs);

			append_objective_with_visual_score(idx_mp, solver);

			//fix s
			std::vector<std::vector<int>> tmpvid = { { 0 } };
			std::vector<std::vector<double>> tmpval = { { -1.0 } };
			solver.append_fixed_var_cons(tmpvid, tmpval);

		}
		else {
			solver.append_vars(n_opt_vars, 0);

			append_orientation_cons_curve_based(idx_mp, solver);

			append_x_regular_cons_curve_based(idx_mp, solver);

			append_imcp_confs(idx_mp, solver, imcompatible_confs);

			//fix s
			std::vector<std::vector<int>> tmpvid = { { 0 } };
			std::vector<std::vector<double>> tmpval = { { -1.0 } };
			solver.append_fixed_var_cons(tmpvid, tmpval);

			append_objective_with_visual_score(idx_mp, solver);
		
		}
	}

	inline std::pair<BlackMesh::BlackMesh<double>*,bool> solve(std::vector<std::vector<std::vector<int>>> &ccid_footprint) {//interation - newccid -old ccid, true-finished
		std::map<int, int> node; std::vector<std::vector<double>> result_solution;
		//solve_problem_with_fixed_solution_capi(node, result_solution);
		std::vector<std::vector<double>> xopt;
		branch_and_bound(xopt);
		if (xopt.size() == 0) {
			std::cout << "Cannot find solution\n";
#if defined(_RECORD_ERROR_LOGS)
			(*_log_off) << _file_under_processing << " ERROR cannot find solution\n";
			(*_log_off).flush();
#endif
			return std::make_pair((BlackMesh::BlackMesh<double>*)NULL,true);
		}

		//show_solution(xopt);

		std::vector<int> final_conf;
		turn_conf_to_stacked(xopt, final_conf);
		write_mesh_change_topo_with_solution(final_conf);
		//write_mesh_grouped_with_solution(final_conf);
		std::vector<bool> del_ls;
		bool all_good=postProcess(final_conf, DELETE_ALL_IN, del_ls, ccid_footprint);

		if (all_good) {
			return std::make_pair(_msh, true);
		}

		return std::make_pair(_msh,false);
	}
	void Optimizer::branch_and_bound(std::vector<std::vector<double>> &xopt);
	double solve_problem_with_fixed_solution(std::map<int, int>& node, std::vector<std::vector<double>>& result_solution);

	
	////build variable map, from (edgeid,confid)->variable id, s will be the idx0 varible with key(-1,-1)
	inline void build_variable_idx_map(map<std::pair<int, int>, int> &mp, const std::map<int, int>& node) {
		mp.clear();
		mp.insert(std::make_pair(std::make_pair(-1,-1), 0));//s

		int counter = 1;
		for (int i = 0; i < _gh.edge_list.size(); i++) {
			////if (node.find(i) != node.end())continue;//fixed node will not be treated as varibles
			//fixed nodes will be added as extra cons
			int degree = _gh.edge_list[i].incident_components.size();
			int siz_conf = configuration_list_all_degrees[degree - 2].size();
			for (int j = 0; j < siz_conf; j++) {
				mp.insert(std::make_pair(std::make_pair(i, j), counter));
				counter++;
			}
		}

		n_opt_vars = counter;
	}
	inline void build_variable_idx_map(map<std::pair<int, int>, int> &mp) {
		mp.clear();
		mp.insert(std::make_pair(std::make_pair(-1, -1), 0));//s

		int counter = 1;
		for (int i = 0; i < _gh.edge_list.size(); i++) {
			////if (node.find(i) != node.end())continue;//fixed node will not be treated as varibles
			//fixed nodes will be added as extra cons
			int degree = _gh.edge_list[i].incident_components.size();
			int siz_conf = configuration_list_all_degrees[degree - 2].size();
			for (int j = 0; j < siz_conf; j++) {
				mp.insert(std::make_pair(std::make_pair(i, j), counter));
				counter++;
			}
		}

		n_opt_vars = counter;
	}
	inline void build_variable_idx_map_curve_based(map<std::pair<int, int>, int> &mp) {
		mp.clear();
		mp.insert(std::make_pair(std::make_pair(-1, -1), 0));//s

		int counter = 1;
		for (int i = 0; i < _gh.non_manifold_curves.size(); i++) {
			////if (node.find(i) != node.end())continue;//fixed node will not be treated as varibles
			//fixed nodes will be added as extra cons

			int degree = _gh.edge_list[_gh.non_manifold_curves[i][0]].incident_components.size();
			int siz_conf = configuration_list_all_degrees[degree - 2].size();
			for (int j = 0; j < siz_conf; j++) {
				mp.insert(std::make_pair(std::make_pair(i, j), counter));
				counter++;
			}
		}

		n_opt_vars = counter;
	}
	////generated mosek wrapper needed data
	//ccid -the lap index of each constraint
	//vid-cons id, variable id
	//val-cons id, coeff

	void Optimizer::append_orientation_cons(map<std::pair<int, int>, int> &mp, Mosek_c_wrapper &solver);
	void Optimizer::append_orientation_cons_curve_based(map<std::pair<int, int>, int> &mp, Mosek_c_wrapper &solver);

	void append_x_regular_cons(map<std::pair<int, int>, int> &mp, Mosek_c_wrapper &solver);
	void Optimizer::append_x_regular_cons_curve_based(map<std::pair<int, int>, int>& mp, Mosek_c_wrapper &solver);

	std::array<double, 5> solve_problem_with_fixed_solution_capi(std::map<int, int>& node, std::vector<std::vector<double>>& result_solution);

	void Optimizer::append_global_cons_with_all_fixed(map<std::pair<int, int>, int> &mp, std::vector<double> solution, Mosek_c_wrapper*);//to compute bound

	void append_fixed_variables(map<std::pair<int, int>, int> &mp, Mosek_c_wrapper & solver, std::map<int, int> node);

	void append_fixed_variables_debug_use(Mosek_c_wrapper & solver, std::vector<double> &s);

	//rounding
	void  Optimizer::maximal_rounding(std::vector<std::vector<double>>& input, std::vector<std::vector<double>>& output);

	//feassible
	bool Optimizer::if_sum_x_equal_to_one(std::vector<std::vector<double>>& x);
	bool Optimizer::if_x_meet_orientation_cons(std::vector<std::vector<double>>& x);
	bool Optimizer::if_x_feasible(std::vector<std::vector<double>>& x);
	bool Optimizer::if_branched_solution_meet_cons(std::map<int, int> node);

	//energy
	double Optimizer::algebraic_conn_bound(std::vector<std::vector<double>>& x);
	std::array<double,5> Optimizer::algebraic_conn_bound_capi(std::vector<std::vector<double>>& x);
	std::array<double, 5> Optimizer::algebraic_conn_bound_capi_debug_use(std::vector<std::vector<double>>& x);
	
	//score
	void update_unary_and_binary_score(double CUTTING_EPS=1e-6);
	void append_objective_with_visual_score(map<std::pair<int, int>, int>& mp, Mosek_c_wrapper &solver);

	//write mesh
	void write_mesh_grouped_with_solution(std::vector<int> &x);
	void write_mesh_change_topo_with_solution(std::vector<int> &x);
	void Optimizer::get_directioend_component_list(std::vector<int> &x, std::vector<int> &ccdir, std::vector<std::vector<int>> &ncclist);

	//change mesh according to solution
	bool Optimizer::postProcess(std::vector<int>& x, POSTPRO_TYPE post_type,
		std::vector<bool> & del_cc_ls, std::vector<std::vector<std::vector<int>>> &ccid_footprint);
	void Optimizer::change_mesh_topo_with_solution(std::vector<int>& x, std::vector<bool> & del_cc_ls);

	//check use, for debug
	bool check_global_laplacian();
	bool check_ori_cons_with_confid();
	void show_solution(std::vector<std::vector<double>> &x);

	//weight update
	void setWeight(double sWeight, double uWeight, double bWeight) {
		_sWeight = sWeight;
		_uWeight = uWeight;
		_bWeight = bWeight;
	};
public:
	Optimizer();
	~Optimizer();
};

