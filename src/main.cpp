#include "stdafx.h"
#include "MeshViewerExt.h"
#include <string>
#include <iostream>
#ifdef _DIR_MODE
#include "dirent.h"
#endif
#include "Utils.h"



template <class Real>
bool batch_preprocess(const char* filename, BlackMesh::BlackMesh<Real> *mesh_, double Z_CUTTING_EPS = 1e-5, double MINIMAL_BORDER = -1) {
	//Z_CUTTING_EPS 1e-5
	//MINIMAL_BORDER -1 means disable
	Preprocess<Real> precesser;
	Preprocess<Real> *precesser_ = &precesser;

	string fm(filename);
	auto found = fm.find_last_of("/\\");
	string filem = fm.substr(found + 1);
	found = filem.find_last_of(".");
	filem = filem.substr(0, found);
	char buffer[512];

	if (_file_under_processing == NULL) {
		_file_under_processing = new char[512];
		sprintf(_file_under_processing, "%s", filem.c_str());
	}

	found = fm.find_last_of(".");
	filem= fm.substr(0, found);

	BlackMesh::BlackMesh<double> shadowMesh;
	shadowMesh.construct_from<ExactScalar>(*mesh_);

	shadowMesh.scale_to_unit();
	mesh_->construct_from<double>(shadowMesh);

	precesser_->set_mesh(mesh_);
	precesser_->opt_vert_ext_on_pro(&shadowMesh,Real(1e-10), false, true);

	precesser_->delete_duplicate_faces();

	shadowMesh.construct_from<ExactScalar>(*mesh_);


	Tribunal tri;
	MeshSampling sampler;

	std::vector<double> posScore, negScore;

	sampler.setMesh(&shadowMesh);


	int cc_deleted = -1;

	for (int i = 0; i < _MAX_DEL_ITER_NUM; i++) {
		sampler.generate_area_and_other_things();
		mesh_->mark_component_with_coherence();
		tri.setMesh(&shadowMesh);

		tri.setSamplig(&sampler);
		tri.setScoreList(NULL, NULL, NULL);
		posScore.resize(sampler.cc_triid.size(), 0);
		negScore.resize(sampler.cc_triid.size(), 0);
		//update orientation score
		tri.update_unary_score(1e-6, 1e-7, &(posScore), &(negScore));

		/*******make coherency******/
		precesser_->make_cc_coherent(posScore, negScore, shadowMesh);


#ifdef _OUTPUT_COHERENT_MESH
		sprintf(buffer, "Result\\%s_debug_coherent_mesh_iter_%d.obj", filem.c_str(), i);
		MeshIO::writeOBJ<Real>(buffer, *(mesh_));
#endif 



		/*******generate sampling points******/
		shadowMesh.clearNumComponent();
		mesh_->clearNumComponent();
		mesh_->mark_component_with_coherence();
		sampler.generate_sampling_points(1000, 20);

		/*******viusal deleting******/
		VisualProcesser vp;
		vp.setMesh(&shadowMesh);

		vp.setSampler(&sampler);
		std::vector<bool> cc_keep_list;

		// offset scale; IMG EPS, zFar, minimal bourder ,CUTTING_PLANE, output

		//	double minmal_boder = iteration_ct < 2 ? -1 : MINIMAL_BORDER;
		double minmal_boder = MINIMAL_BORDER;

#ifdef _MAX_CC_NM_FOR_PREPROCESSING
		if (mesh_->GetNumComponents() > _MAX_CC_NM_FOR_PREPROCESSING) {
#if defined(_RECORD_ERROR_LOGS)
			(*_log_off) << _file_under_processing << "ERROR too many ccs in preprocessing\n";
			(*_log_off).flush();
#endif
			return false;
		}
#endif

		cc_deleted = vp.delete_duplicate_cc_visual_based(1.0, _eps_vis, 2.8f, minmal_boder, Z_CUTTING_EPS, cc_keep_list);
		
#ifdef _OUTPUT_VD_DEL_AS_MESH
		vp.save_res_mesh(cc_keep_list, i);
#endif

		mesh_->update_mesh_to_keep_list(cc_keep_list);
		shadowMesh.update_mesh_to_keep_list(cc_keep_list);


		if (i == 0) {

			/******* Intersection, merge, del duplicate******/


#if defined(_RECORD_ERROR_LOGS)
			if(mesh_->GetNumEdges()!=shadowMesh.GetNumEdges()|| 
				mesh_->GetNumVertices() != shadowMesh.GetNumVertices()||
				mesh_->GetNumTriangles() != shadowMesh.GetNumTriangles()){

			(*_log_off) << _file_under_processing << "ERROR shadowMesh not match\n";
			(*_log_off).flush();
			}
#endif


			Intersector<double, Real> intsc;
			intsc.run(&shadowMesh, mesh_);

			//post process
			precesser_->set_mesh(mesh_);
			precesser_->opt_vert_ext(Real(0), true);

			precesser_->delete_duplicate_faces();

			shadowMesh.construct_from<ExactScalar>(*mesh_);


#ifdef _OUTPUT_RESOLVED_MESH
			sprintf(buffer, "Result\\%s_debug_resolved_mesh_iter%d_merged.obj", filem.c_str(), i);
			MeshIO::writeOBJ<Real>(buffer, *(mesh_));
#endif
		}
	}

	if (cc_deleted > 0) {

		sampler.generate_area_and_other_things();
		mesh_->mark_component_with_coherence();
		tri.setMesh(&shadowMesh);

		tri.setSamplig(&sampler);
		tri.setScoreList(NULL, NULL, NULL);
		posScore.resize(sampler.cc_triid.size(), 0);
		negScore.resize(sampler.cc_triid.size(), 0);
		//update orientation score
		tri.update_unary_score(1e-6, 1e-7, &(posScore), &(negScore));

		/*******make coherency******/
		precesser_->make_cc_coherent(posScore, negScore, shadowMesh);

	}
#ifdef _DIR_MODE
	sprintf(buffer, "Result\\%s_vd.obj", filem.c_str());
#else
	sprintf(buffer, "%s_vd.obj", filem.c_str());
#endif
	MeshIO::writeOBJ<Real>(buffer, *(mesh_));


	return true;
}




void batch_process(const char * filename, BlackMesh::BlackMesh<double> **msh) {


			Optimizer opt;

			//for (int i = 0; i < res.size(); i++) {
			opt.setMesh(*msh);
			opt.setTimer(_tmer, 1000);
			opt.setWeight(_sw, _uw, _bw);


			string fm(filename);

			auto found = fm.find_last_of("/\\");
			string filem = fm.substr(found + 1);
			found = filem.find_last_of(".");
			filem = filem.substr(0, found);

			if (_file_under_processing == NULL) {
				_file_under_processing = new char[512];
				sprintf(_file_under_processing, "%s", filem.c_str());
			}

#ifndef _DIR_MODE
			 //found = fm.find_last_of(".");
			// filem = fm.substr(0, found);
			found= fm.find_last_of("/\\");
			string pathnm = fm.substr(0, found);
#endif

			string laptype;
			string bintype;

			laptype = "cvlp";
			bintype = "cvbin";


			std::vector<std::vector<std::vector<int>>> cc_ftpt;
			for (int iterct = 0; iterct < 2; iterct++) {

				int succ = opt.pre_assemble_curve_based();

				if (succ == 0) {


					char buffer[128];
					if (iterct == 0) {
						sprintf(buffer, "%s_opt", filem.c_str());
					}
					else {
						sprintf(buffer, "%s_opt_HR", filem.c_str());
					}

					//if (iterct == 0) {
					//	char filemaneNoModi[128];
					//	sprintf(filemaneNoModi, "Result\\%sO.obj", buffer);
					//	MeshIO::writeOBJ<double>(filemaneNoModi, **msh);
					//}
#ifdef _DIR_MODE
					opt.setOutputName("Result\\", buffer);
#else
					char pth[512];
					sprintf(pth, "%s", pathnm.c_str());
					opt.setOutputName(pth, buffer);
#endif
					std::cout << "Solving Mesh\n";
					auto ptr = opt.solve(cc_ftpt);
					if (ptr.first != NULL) {
						(*msh) = ptr.first;
					}

					if (ptr.second) {
						break;
					}
			}
				else {

					char buffer[512];
					if (iterct == 0) {
						sprintf(buffer, "%s_opt", filem.c_str());
					}
					else {
						sprintf(buffer, "%s_opt_HR", filem.c_str());
					}

					char filemaneNoModi[512];
#ifdef _DIR_MODE
					sprintf(filemaneNoModi, "Result\\%sN.obj", buffer);
#else
					sprintf(filemaneNoModi, "%s%sN.obj", pathnm.c_str(), buffer);
#endif
					MeshIO::writeOBJ<double>(filemaneNoModi, **(msh));
					std::cout << "Mesh does not need opt\n";
					break;
				}
			}


	std::cout << "Task Finished\n";
}


void run(const char* filename, double Z_CUTTING_EPS, double minimal_border) {
#ifdef _DIR_MODE
	{// extract folder, set as current path
		string tmp(filename);
		auto found = tmp.find_last_of("/\\");
		string tpath = tmp.substr(0, found);

		int  size = tpath.size();
		wchar_t *bufferff = new wchar_t[size + 1];
		MultiByteToWideChar(CP_ACP, 0, tpath.c_str(), size, bufferff, size * sizeof(wchar_t));
		bufferff[size] = 0;

		SetCurrentDirectory(bufferff);
	}
#endif
	BlackMesh::BlackMesh<ExactScalar> *msh;
	msh = new BlackMesh::BlackMesh<ExactScalar>();
	bool succ = MeshIO::readMesh(filename, *msh);


	if (succ = false)
	{
		std::cout << "cannot read";
		return;
	}

	bool continue_flag = true;
#ifndef _DISABLE_VD
	continue_flag = batch_preprocess(filename, msh, Z_CUTTING_EPS, minimal_border);
#endif//_DISABLE_VD
#ifndef _DISABLE_OPT
	if (continue_flag) {
#ifdef _DISABLE_VD
		msh->scale_to_unit();
#endif//_DISABLE_VD

		BlackMesh::BlackMesh<double> *shadowMesh;
		shadowMesh = new BlackMesh::BlackMesh<double>();
		shadowMesh->construct_from<ExactScalar>(*msh);
		batch_process(filename, &shadowMesh);

		delete shadowMesh;

}
#endif//_DISABLE_OPT

	delete msh;
}


inline bool exists_fils(const std::string& name) {
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}

int main(int argc, char* argv[]) {
	if (argc == 7) {

		//filename eps_vis cslmin w1 w2 #threads

		glutInit(&argc, argv);
		MeshViewerExt viewer("Mesh Viewer", 1000, 600);
		viewer.launch();
		//glutHideWindow();

		double zcutting = 1e-5;
		_eps_vis = atof(argv[2]);
		double border = atof(argv[3]);

		_sw = -1;
		_uw = atof(argv[4]);
		_bw = atof(argv[5]);

		_opt_max_iter = 1;

		_tmer = 5;

		_num_threads = atoi(argv[6]);

		omp_set_num_threads(_num_threads);

		printf("Eps_vis: %.8f\nCslMin: %.8f\nw1: %.8f\nw2: %.8f\n#threads: %d\n",_eps_vis,border, _uw, _bw,_num_threads);

#if defined(_RECORD_ERROR_LOGS) 
		_log_off = new ofstream("_DEBUG_LOGS_ALL.txt", ios::app);
#endif

#ifdef _DIR_MODE

		DIR *dir;
		struct dirent *ent;
		if ((dir = opendir(argv[1])) != NULL) {
			/* print all the files and directories within directory */
			string folderpath(argv[1]);
			while ((ent = readdir(dir)) != NULL) {
				printf("@@@Processing %s\n", ent->d_name);
				std::string filename(ent->d_name);
				auto found = filename.find_last_of(".");
				string filesuffix = filename.substr(found + 1);
				if (filesuffix.compare("off") == 0 || filesuffix.compare("obj") == 0) {

					//check if already run
					
					string testfileWholePath = folderpath;
					string filewosuffix= filename.substr(0,found);
					testfileWholePath.append("Result\\");
					testfileWholePath.append(filewosuffix);
#ifndef _DISABLE_VD
					testfileWholePath.append("_vd.obj");
#else
					testfileWholePath.append("_opt.obj");
#endif

					if(!exists_fils(testfileWholePath)){

					string fileWholePath = folderpath;
					fileWholePath.append(filename);


					try {
						srand(0);
						run(fileWholePath.c_str(), zcutting, border);
					}
					catch (...) {
						(*_log_off) << filename << " Error, Crashed\n";
					
					}
					delete[]_file_under_processing;
					_file_under_processing = NULL;
					}
					else {
						std::cout << "**Already Processed\n";
					}
				}
			}
			closedir(dir);
		}
		else {
			/* could not open directory */
			perror("");
			return EXIT_FAILURE;
		}

#else
		printf("@@@Processing %s\n", argv[1]);
		std::string filename(argv[1]);
		auto found = filename.find_last_of(".");
		string filesuffix = filename.substr(found + 1);

		found = filename.find_last_of("\\");
		string folderpath = filename.substr(0,found);
		if (filesuffix.compare("off") == 0 || filesuffix.compare("obj") == 0) {


				try {
					srand(0);
					//srand(time(NULL));
					run(filename.c_str(), zcutting, border);
				}
				catch (...) {
					(*_log_off) << filename << " Error, Crashed\n";

				}
				delete[]_file_under_processing;
				_file_under_processing = NULL;
\
	}

#endif

#if defined(_RECORD_ERROR_LOGS)

		delete _log_off;
		_log_off = NULL;
#endif

		//system("Pause");
}
	else {
		std::cout << "Not enough para\n";
	}

	return 0;
}