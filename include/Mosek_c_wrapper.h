#pragma once
#include "mosek.h"
#include <vector>
#include <iostream>
class Mosek_c_wrapper
{
private:
	MSKrescodee  r;
	MSKenv_t     env = NULL;
	MSKtask_t    task = NULL;
	double falpha = 1.0;
	
public:
	Mosek_c_wrapper() {
		//r = MSK_makeenv(&env, NULL);
		//make_task;
	}

	inline MSKenv_t get_Env() { return env; }
	inline MSKtask_t getTask() { return task; }

	inline bool make_task() {

		//MSKrescodee(MSKAPI MSK_maketask) (
		//	MSKenv_t env,
		//	MSKint32t maxnumcon,
		//	MSKint32t maxnumvar,
		//	MSKtask_t * task)
		//	Creates a new task.

		//	Parameters:

		//  env(MSKenv_t) ?The MOSEK environment. (input)
		//	maxnumcon(MSKint32t) ?An optional estimate on the maximum number of constraints in the task.Can be 0 if no such estimate is known. (input)
		//	maxnumvar(MSKint32t) ?An optional estimate on the maximum number of variables in the task.Can be 0 if no such estimate is known. (input)
		//	task(MSKtask_t by reference) ?An optimization task. (output)
		r = MSK_makeenv(&env, NULL);
		if (r != MSK_RES_OK) return false;
		r = MSK_maketask(env, 0, 0, &task);
#ifdef _OUTPUT_MOSEK_LOG
		if (r == MSK_RES_OK)
		r = MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);
#endif
#ifdef _MAX_TIME_FOR_MOSEK
		setMaxOptSolverTime(_MAX_TIME_FOR_MOSEK);
#endif
		return true;
	}

	inline void setMaxOptSolverTime(double tm) {//mins
		if (r == MSK_RES_OK)
		r= MSK_putnadouparam(task, "MSK_DPAR_OPTIMIZER_MAX_TIME", tm*60.0);
	}

	inline int get_n_vars() {
		int varidx;
		if (r == MSK_RES_OK)
			r = MSK_getnumvar(task, &varidx);
		return varidx;
	}

	inline void append_vars(int NVARS, int SDP_dim) {
		//one unbounded and others from [0,1], one SDP
		if (r == MSK_RES_OK)
		{
			r = MSK_appendvars(task, NVARS);

			//s
			if (r == MSK_RES_OK) {
				r = MSK_putvarbound(task,
					0,
					MSK_BK_FR,
					-MSK_INFINITY,
					MSK_INFINITY);
			}
			
			//confs
			for (MSKint32t j = 1; j < NVARS && r == MSK_RES_OK; ++j) {
				r = MSK_putvarbound(task,
					j,
					MSK_BK_RA,
					0.0,
					1.0);
			}
		}

		if(SDP_dim>0)
		append_one_SDP_var(SDP_dim);
	}

	inline void append_one_SDP_var(int dim) {
		if (r == MSK_RES_OK) {
			r = MSK_appendbarvars(task, 1, &MSKint32t(dim));
		}
	}

	inline void append_empty_cons(int NCONS) {
		if (r == MSK_RES_OK)
		{
			r = MSK_appendcons(task, NCONS);
		}
	}

	inline void append_global_connectivity_eq_cons(std::vector<std::pair<int,int>> &ccid, 
		int dim, std::vector<std::vector<int>> &vid, std::vector<std::vector<double>> &val) {
		MSKint32t conidx;
		if (r == MSK_RES_OK) {
			r = MSK_getnumcon(task, &conidx);
			if (conidx != 0) { 
				std::cout << "Initial cons is not zero\n";
#if defined(_RECORD_ERROR_LOGS)
				(*_log_off) << _file_under_processing << " ERROR Initial cons is not zero\n";
				(*_log_off).flush();
#endif
			}
		}

		/* Append  new constraints */
		if (r == MSK_RES_OK)
		{
			r = MSK_appendcons(task, ccid.size());
			//numcon++;
		}

		/* Set the bounds on constraints.*/
		int current_con_siz;
		if (r == MSK_RES_OK)
			r = MSK_getnumcon(task, &current_con_siz);
		for (int i = conidx; i < current_con_siz && r == MSK_RES_OK; ++i)
			r = MSK_putconbound(task,
				i,
				MSK_BK_FX,
				0.0,
				0.0);

		//set a
		for (int i = conidx; i < current_con_siz && r == MSK_RES_OK; ++i) {
			/*****
			MSKrescodee(MSKAPI MSK_putarow) (
				MSKtask_t task,
				MSKint32t i,
				MSKint32t nzi,
				const MSKint32t * subi,
				const MSKrealt * vali)

				task (MSKtask_t) ?An optimization task. (input)
				i (MSKint32t) ?Index of a row in A. (input)
				nzi (MSKint32t) ?Number of non-zeros in row i of A. (input)
				subi (MSKint32t*) ?Column indexes of non-zero values in row i of A. (input)
				vali (MSKrealt*) ?New non-zero values of row i in A. (input)
				***/
			//need to make sure the vector wont be released during using!!!
			r = MSK_putarow(task,
				i,
				vid[i].size(),
				vid[i].data(),
				val[i].data());
		}

		/* Add the rows of barA */
		for (int i = conidx; i < current_con_siz && r == MSK_RES_OK; ++i) {
		/*****
		MSKrescodee (MSKAPI MSK_appendsparsesymmat) (
			MSKtask_t task,
			MSKint32t dim,
			MSKint64t nz,
			const MSKint32t * subi,
			const MSKint32t * subj,
			const MSKrealt * valij,
			MSKint64t * idx)

			task (MSKtask_t) ?An optimization task. (input)
			dim (MSKint32t) ?Dimension of the symmetric matrix that is appended. (input)
			nz (MSKint64t) ?Number of triplets. (input)
			subi (MSKint32t*) ?Row subscript in the triplets. (input)
			subj (MSKint32t*) ?Column subscripts in the triplets. (input)
			valij (MSKrealt*) ?Values of each triplet. (input)
			idx (MSKint64t by reference) ?Unique index assigned to the inputted matrix that can be used for later reference. (output)
		*****/
			MSKint64t idx;
			MSKint32t rid[] = { ccid[i].second };
			MSKint32t cid[] = { ccid[i].first };
			MSKrealt vv[] = { -1.0 };
			r = MSK_appendsparsesymmat(task,
				dim,
				1,
				rid,
				cid,
				vv,
				&idx);
			/************
			MSKrescodee (MSKAPI MSK_putbaraij) (
					MSKtask_t task,
					MSKint32t i,
					MSKint32t j,
					MSKint64t num,
					const MSKint64t * sub,
					const MSKrealt * weights)

					task (MSKtask_t) ?An optimization task. (input)
					i (MSKint32t) ?Row index of A¯¯¯¯. (input)
					j (MSKint32t) ?Column index of A¯¯¯¯. (input)
					num (MSKint64t) ?The number of terms in the weighted sum that forms A¯¯¯¯ij. (input)???????????????
					sub (MSKint64t*) ?Indices in E of the matrices appearing in the weighted sum for A¯¯¯¯ij. (input)???????????????
					weights (MSKrealt*) ?weights[k] is the coefficient of the sub[k]-th element of E in the weighted sum forming A¯¯¯¯ij. (input)
			**********/

			if (r == MSK_RES_OK)
				r = MSK_putbaraij(task, i, 0, 1, &idx, &falpha);
		}
	}

	inline void append_orientation_cons(std::vector<std::vector<int>> &vid, std::vector<std::vector<double>> &val) {
		MSKint32t conidx;
		if (r == MSK_RES_OK) {
			r = MSK_getnumcon(task, &conidx);
			//if (conidx != 0)std::cout << "Initial cons is not zero\n";
		}

		/* Append  new constraints */
		if (r == MSK_RES_OK)
		{
			r = MSK_appendcons(task, vid.size());
			//numcon++;
		}

		/* Set the bounds on constraints.*/
		int current_con_siz;
		if (r == MSK_RES_OK)
			r = MSK_getnumcon(task, &current_con_siz);
		for (int i = conidx; i < current_con_siz && r == MSK_RES_OK; ++i)
			r = MSK_putconbound(task,
				i,
				MSK_BK_FX,
				0.0,
				0.0);

		//set a
		int rowid = 0;
		for (int i = conidx; i < current_con_siz && r == MSK_RES_OK; ++i) {
			/*****
			MSKrescodee(MSKAPI MSK_putarow) (
			MSKtask_t task,
			MSKint32t i,
			MSKint32t nzi,
			const MSKint32t * subi,
			const MSKrealt * vali)

			task (MSKtask_t) ?An optimization task. (input)
			i (MSKint32t) ?Index of a row in A. (input)
			nzi (MSKint32t) ?Number of non-zeros in row i of A. (input)
			subi (MSKint32t*) ?Column indexes of non-zero values in row i of A. (input)
			vali (MSKrealt*) ?New non-zero values of row i in A. (input)
			***/
			//need to make sure the vector wont be released during using!!!
			r = MSK_putarow(task,
				i,
				vid[rowid].size(),
				vid[rowid].data(),
				val[rowid].data());
			rowid++;
		}
	}


	inline void append_regular_cons(std::vector<std::vector<int>> &vid, std::vector<std::vector<double>> &val) {
		MSKint32t conidx;
		if (r == MSK_RES_OK) {
			r = MSK_getnumcon(task, &conidx);
			//if (conidx != 0)std::cout << "Initial cons is not zero\n";
		}

		/* Append  new constraints */
		if (r == MSK_RES_OK)
		{
			r = MSK_appendcons(task, vid.size());
			//numcon++;
		}

		/* Set the bounds on constraints.*/
		int current_con_siz;
		if (r == MSK_RES_OK)
			r = MSK_getnumcon(task, &current_con_siz);
		for (int i = conidx; i < current_con_siz && r == MSK_RES_OK; ++i)
			r = MSK_putconbound(task,
				i,
				MSK_BK_FX,
				1.0,
				1.0);

		//set a
		int rowid = 0;
		for (int i = conidx; i < current_con_siz && r == MSK_RES_OK; ++i) {
			/*****
			MSKrescodee(MSKAPI MSK_putarow) (
			MSKtask_t task,
			MSKint32t i,
			MSKint32t nzi,
			const MSKint32t * subi,
			const MSKrealt * vali)

			task (MSKtask_t) ?An optimization task. (input)
			i (MSKint32t) ?Index of a row in A. (input)
			nzi (MSKint32t) ?Number of non-zeros in row i of A. (input)
			subi (MSKint32t*) ?Column indexes of non-zero values in row i of A. (input)
			vali (MSKrealt*) ?New non-zero values of row i in A. (input)
			***/
			//need to make sure the vector wont be released during using!!!
			r = MSK_putarow(task,
				i,
				vid[rowid].size(),
				vid[rowid].data(),
				val[rowid].data());
			rowid++;
		}
	}

	inline void append_fixed_var_cons(std::vector<std::vector<int>> &vid, std::vector<std::vector<double>> &val, double rhsterm=1.0) {
		MSKint32t conidx;
		if (r == MSK_RES_OK) {
			r = MSK_getnumcon(task, &conidx);
			//if (conidx != 0)std::cout << "Initial cons is not zero\n";
		}

		/* Append  new constraints */
		if (r == MSK_RES_OK)
		{
			r = MSK_appendcons(task, vid.size());
			//numcon++;
		}

		/* Set the bounds on constraints.*/
		int current_con_siz;
		if (r == MSK_RES_OK)
			r = MSK_getnumcon(task, &current_con_siz);
		for (int i = conidx; i < current_con_siz && r == MSK_RES_OK; ++i)
			r = MSK_putconbound(task,
				i,
				MSK_BK_FX,
				rhsterm,
				rhsterm);

		//set a
		int rowid = 0;
		for (int i = conidx; i < current_con_siz && r == MSK_RES_OK; ++i) {
			/*****
			MSKrescodee(MSKAPI MSK_putarow) (
			MSKtask_t task,
			MSKint32t i,
			MSKint32t nzi,
			const MSKint32t * subi,
			const MSKrealt * vali)

			task (MSKtask_t) ?An optimization task. (input)
			i (MSKint32t) ?Index of a row in A. (input)
			nzi (MSKint32t) ?Number of non-zeros in row i of A. (input)
			subi (MSKint32t*) ?Column indexes of non-zero values in row i of A. (input)
			vali (MSKrealt*) ?New non-zero values of row i in A. (input)
			***/
			//need to make sure the vector wont be released during using!!!
			r = MSK_putarow(task,
				i,
				vid[rowid].size(),
				vid[rowid].data(),
				val[rowid].data());
			rowid++;
		}
	}

	inline void append_fixed_var_cons_debug_use(std::vector<std::vector<int>> &vid, std::vector<std::vector<double>> &val, std::vector<double> &RHS) {
		MSKint32t conidx;
		if (r == MSK_RES_OK) {
			r = MSK_getnumcon(task, &conidx);
			//if (conidx != 0)std::cout << "Initial cons is not zero\n";
		}

		/* Append  new constraints */
		if (r == MSK_RES_OK)
		{
			r = MSK_appendcons(task, vid.size());
			//numcon++;
		}

		/* Set the bounds on constraints.*/
		int current_con_siz;
		if (r == MSK_RES_OK)
			r = MSK_getnumcon(task, &current_con_siz);
		for (int i = conidx; i < current_con_siz && r == MSK_RES_OK; ++i)
			r = MSK_putconbound(task,
				i,
				MSK_BK_FX,
				RHS[i- conidx],
				RHS[i - conidx]);

		//set a
		int rowid = 0;
		for (int i = conidx; i < current_con_siz && r == MSK_RES_OK; ++i) {
			/*****
			MSKrescodee(MSKAPI MSK_putarow) (
			MSKtask_t task,
			MSKint32t i,
			MSKint32t nzi,
			const MSKint32t * subi,
			const MSKrealt * vali)

			task (MSKtask_t) ?An optimization task. (input)
			i (MSKint32t) ?Index of a row in A. (input)
			nzi (MSKint32t) ?Number of non-zeros in row i of A. (input)
			subi (MSKint32t*) ?Column indexes of non-zero values in row i of A. (input)
			vali (MSKrealt*) ?New non-zero values of row i in A. (input)
			***/
			//need to make sure the vector wont be released during using!!!
			r = MSK_putarow(task,
				i,
				vid[rowid].size(),
				vid[rowid].data(),
				val[rowid].data());
			rowid++;
		}
	}

	inline void append_all_fixed_conn_cons(std::vector<std::pair<int, int>> &ccid,
		int dim, std::vector<std::vector<int>> &vid, std::vector<std::vector<double>> &val) {
		MSKint32t conidx;
		if (r == MSK_RES_OK) {
			r = MSK_getnumcon(task, &conidx);
			if (conidx != 0) { 
				std::cout << "Initial cons is not zero\n"; 
#if defined(_RECORD_ERROR_LOGS)
				(*_log_off) << _file_under_processing << " ERROR Initial cons is not zero\n";
				(*_log_off).flush();
#endif
			}
		}

		/* Append  new constraints */
		if (r == MSK_RES_OK)
		{
			r = MSK_appendcons(task, ccid.size());
			//numcon++;
		}

		/* Set the bounds on constraints.*/
		int current_con_siz;
		if (r == MSK_RES_OK)
			r = MSK_getnumcon(task, &current_con_siz);
		for (int i = conidx; i < current_con_siz && r == MSK_RES_OK; ++i)
			r = MSK_putconbound(task,
				i,
				MSK_BK_FX,
				0.0,
				0.0);

		//set a
		for (int i = conidx; i < current_con_siz && r == MSK_RES_OK; ++i) {
			/*****
			MSKrescodee(MSKAPI MSK_putarow) (
			MSKtask_t task,
			MSKint32t i,
			MSKint32t nzi,
			const MSKint32t * subi,
			const MSKrealt * vali)

			task (MSKtask_t) ?An optimization task. (input)
			i (MSKint32t) ?Index of a row in A. (input)
			nzi (MSKint32t) ?Number of non-zeros in row i of A. (input)
			subi (MSKint32t*) ?Column indexes of non-zero values in row i of A. (input)
			vali (MSKrealt*) ?New non-zero values of row i in A. (input)
			***/
			//need to make sure the vector wont be released during using!!!
			r = MSK_putarow(task,
				i,
				vid[i].size(),
				vid[i].data(),
				val[i].data());
		}

		/* Add the rows of barA */
		for (int i = conidx; i < current_con_siz && r == MSK_RES_OK; ++i) {
			/*****
			MSKrescodee (MSKAPI MSK_appendsparsesymmat) (
			MSKtask_t task,
			MSKint32t dim,
			MSKint64t nz,
			const MSKint32t * subi,
			const MSKint32t * subj,
			const MSKrealt * valij,
			MSKint64t * idx)

			task (MSKtask_t) ?An optimization task. (input)
			dim (MSKint32t) ?Dimension of the symmetric matrix that is appended. (input)
			nz (MSKint64t) ?Number of triplets. (input)
			subi (MSKint32t*) ?Row subscript in the triplets. (input)
			subj (MSKint32t*) ?Column subscripts in the triplets. (input)
			valij (MSKrealt*) ?Values of each triplet. (input)
			idx (MSKint64t by reference) ?Unique index assigned to the inputted matrix that can be used for later reference. (output)
			*****/
			MSKint64t idx;
			MSKint32t rid[] = { ccid[i].second };
			MSKint32t cid[] = { ccid[i].first };
			MSKrealt vv[] = { -1.0 };
			r = MSK_appendsparsesymmat(task,
				dim,
				1,
				rid,
				cid,
				vv,
				&idx);
			/************
			MSKrescodee (MSKAPI MSK_putbaraij) (
			MSKtask_t task,
			MSKint32t i,
			MSKint32t j,
			MSKint64t num,
			const MSKint64t * sub,
			const MSKrealt * weights)

			task (MSKtask_t) ?An optimization task. (input)
			i (MSKint32t) ?Row index of A¯¯¯¯. (input)
			j (MSKint32t) ?Column index of A¯¯¯¯. (input)
			num (MSKint64t) ?The number of terms in the weighted sum that forms A¯¯¯¯ij. (input)???????????????
			sub (MSKint64t*) ?Indices in E of the matrices appearing in the weighted sum for A¯¯¯¯ij. (input)???????????????
			weights (MSKrealt*) ?weights[k] is the coefficient of the sub[k]-th element of E in the weighted sum forming A¯¯¯¯ij. (input)
			**********/

			if (r == MSK_RES_OK)
				r = MSK_putbaraij(task, i, 0, 1, &idx, &falpha);
		}
	}

	inline void append_objective() {
		if (r == MSK_RES_OK)
			r = MSK_putcj(task, 0, 1.0);//min, so coeff of s is -1.0

		if (r == MSK_RES_OK)
			r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE);
	}

	inline void append_objective_with_visual_score(std::vector<int> &vid, std::vector<double> &val, double connect_bound_weight) {
		if (r == MSK_RES_OK)
			r = MSK_putcj(task, 0, connect_bound_weight);

		for (int ss = 0; ss < vid.size(); ss++) {
			if (r == MSK_RES_OK)
				r = MSK_putcj(task, vid[ss], val[ss]);
		}

		if (r == MSK_RES_OK)
			r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE);
	}

	inline bool solve() {
	
		bool succ=false;
		if (r == MSK_RES_OK)
		{
			MSKrescodee trmcode;

			/* Run optimizer */
			r = MSK_optimizetrm(task, &trmcode);

			/* Print a summary containing information
			about the solution for debugging purposes*/
			MSK_solutionsummary(task, MSK_STREAM_MSG);

			if (r == MSK_RES_OK)
			{
				MSKsolstae solsta;

				MSK_getsolsta(task, MSK_SOL_ITR, &solsta);

				switch (solsta)
				{
				case MSK_SOL_STA_OPTIMAL:
				case MSK_SOL_STA_NEAR_OPTIMAL:

					succ = true;
					break;

				case MSK_SOL_STA_DUAL_INFEAS_CER:
				case MSK_SOL_STA_PRIM_INFEAS_CER:
				case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
				case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
#ifdef _OUTPUT_MOSEK_LOG
					printf("Primal or dual infeasibility certificate found.\n");
#endif
#if defined(_RECORD_ERROR_LOGS) 
					(*_log_off) << _file_under_processing << " ERROR Primal or dual infeasibility certificate found\n";
					(*_log_off).flush();
#endif
					succ = false;
					break;

				case MSK_SOL_STA_UNKNOWN:

					printf("The status of the solution could not be determined.\n");
#if defined(_RECORD_ERROR_LOGS)
					(*_log_off) << _file_under_processing << " ERROR The status of the solution could not be determined\n";
					(*_log_off).flush();
#endif
					succ = false;
					break;

				default:

					printf("Other solution status.");
#if defined(_RECORD_ERROR_LOGS) 
					(*_log_off) << _file_under_processing << " ERROR Other solution status\n";
					(*_log_off).flush();
#endif
					succ = false;
					break;
				}
			}
			else
			{
				printf("Error while optimizing.\n");
#if defined(_RECORD_ERROR_LOGS)
				(*_log_off) << _file_under_processing << " ERROR Error while optimizing\n";
				(*_log_off).flush();
#endif
				succ = false;
			}
		}

		if (r != MSK_RES_OK)
		{
			/* In case of an error print error code and description. */
			char symname[MSK_MAX_STR_LEN];
			char desc[MSK_MAX_STR_LEN];

			printf("An error occurred while optimizing.\n");
#if defined(_RECORD_ERROR_LOGS)
			(*_log_off) << _file_under_processing << " ERROR An error occurred while optimizing\n";
			(*_log_off).flush();
#endif
			MSK_getcodedesc(r,
				symname,
				desc);
			printf("Error %s - '%s'\n", symname, desc);
			succ = false;
		}
	
		return succ;
	}

	inline void get_solution(double &s, std::vector<double> &x, int NUMVAR) {
		x.clear();
		MSKrealt     *xx;
		xx = (MSKrealt *)MSK_calloctask(task, NUMVAR, sizeof(MSKrealt));
		MSK_getxx(task,
			MSK_SOL_ITR,
			xx);
		s = xx[0];
		for (int i = 1; i < NUMVAR; i++) {
			x.push_back(xx[i]);
		}
		//free(xx);
	}

	inline void get_solution_matrix(std::vector<double> &x, int NUMVAR) {
		x.clear();
		int NUM_MTX = NUMVAR*(NUMVAR + 1) / 2;
		MSKrealt     *barx;
		barx = (MSKrealt *)MSK_calloctask(task, NUM_MTX, sizeof(MSKrealt));

		MSK_getbarxj(task,
			MSK_SOL_ITR,    /* Request the interior solution. */
			0,
			barx);

		for (int i = 0; i < NUM_MTX; i++) {
			x.push_back(barx[i]);
		}

	//	free(barx);
	}

	inline void cloned_from(Mosek_c_wrapper& input) {
		r=MSK_clonetask(input.getTask(), &task);
		env=input.get_Env();
#ifdef _OUTPUT_MOSEK_LOG
		if (r == MSK_RES_OK)
			r = MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);
#endif
	}

	static void MSKAPI printstr(void *handle,
		const char str[])
	{
		printf("%s", str);
	}

	~Mosek_c_wrapper() {
		if (r == MSK_RES_OK)
		MSK_deletetask(&task);
		//MSK_deleteenv(&env);
	}
};

