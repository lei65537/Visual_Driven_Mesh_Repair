#include "stdafx.h"
#include "VisualProcesser.h"

#define _BORDER_TO_BBX 1.2



void VisualProcesser::compute_the_diff_non_blank(std::vector<float>& img1, std::vector<float>& img2,
	std::vector<int>& IDmask, int ccid,
	const  float& clear_val, float& acc_err, int& acc_msk_siz, int image_wid) {//image 2 is the one with cc


																//float err = 0;
	std::vector<float> acc_err_row;
	std::vector<int> acc_msk_siz_row;
	acc_err_row.resize(image_wid, 0.0f);
	acc_msk_siz_row.resize(image_wid, 0);


#pragma omp parallel for
	for (int i = 0; i < image_wid; i++) {
		for (int j = 0; j < image_wid; j++) {
			int eid = i* image_wid + j;
			if (IDmask[eid] == 0)continue;
			if ((IDmask[eid] - 1) != ccid)continue;
			acc_err_row[i] += abs(img1[eid] - img2[eid]);
			acc_msk_siz_row[i]++;
		}
	}

	for (int i = 0; i < image_wid; i++) {
		acc_err += acc_err_row[i];
		acc_msk_siz += acc_msk_siz_row[i];
	}




}


void VisualProcesser::update_score_and_relation(OffScreenRenderer &rdr, float offset_scale, float zFar, float minimal_border,
	float CutEPS, std::set<unsigned int> &update_list, std::set<unsigned int> &delete_list,
	std::vector<std::pair<unsigned int, float>>& score_list, std::vector<std::set<unsigned int>> &relation_list,
	std::vector<std::vector<int>> &IDmask, std::vector<int> &temp_buffer,
	std::vector<std::vector<float>> &image_stack_wt, std::vector<std::vector<float>> &image_stack_wo
) {
	
	int img_wid = rdr.getWidth();

	if (update_list.size() == 0) {
		//the first iteration//all show
		score_list.clear();
		score_list.reserve(msh->GetNumComponents());

		relation_list.clear();
		relation_list.resize(msh->GetNumComponents(), std::set<unsigned int>{});

		for (unsigned int i = 0; i <  sampler->cc_sampled_triid.size(); i++) {
			//i is ccid
			//i = 537;
			//rdr.make_cc_NOTshow(std::set<unsigned int>{536});

			rdr.make_cc_NOTshow_for_compare(i);

			float ccborder = get_cc_border(sampler->cc_bbx[i]);
			ccborder *= _BORDER_TO_BBX;

			if (minimal_border>0&&ccborder < minimal_border)ccborder = minimal_border;

			float zNear = 0.1;

#ifdef _SAVE_OFFRENDER_TO_IMG
			rdr.output_counter = i;
			rdr.rder_couter = 0;
#endif


			std::vector<bool> cc_in_current_scene;
			cc_in_current_scene.resize(sampler->cc_sampled_triid.size(), false);
			float score_current_cc = 0.0f;
			int total_pix_siz = 0;

			for (int j = 0; j <  sampler->cc_sampled_triid[i].size(); j++) {
				//for each triangle in cci

				//assemble render para
				RenderStrc render_para;

				int triid = sampler->cc_sampled_triid[i][j].first;
				auto nm = msh->GetTriangles()[triid].normal;

				render_para.campos_pos[0] = offset_scale*nm[0] + sampler->cc_sampled_triid[i][j].second[0];
				render_para.campos_pos[1] = offset_scale*nm[1] + sampler->cc_sampled_triid[i][j].second[1];
				render_para.campos_pos[2] = offset_scale*nm[2] + sampler->cc_sampled_triid[i][j].second[2];

				render_para.campos_neg[0] = -offset_scale*nm[0] + sampler->cc_sampled_triid[i][j].second[0];
				render_para.campos_neg[1] = -offset_scale*nm[1] + sampler->cc_sampled_triid[i][j].second[1];
				render_para.campos_neg[2] = -offset_scale*nm[2] + sampler->cc_sampled_triid[i][j].second[2];

				render_para.target[0] = sampler->cc_sampled_triid[i][j].second[0];
				render_para.target[1] = sampler->cc_sampled_triid[i][j].second[1];
				render_para.target[2] = sampler->cc_sampled_triid[i][j].second[2];

				std::vector<float> tmp(3);
				tmp[0] = nm[0];
				tmp[1] = nm[1];
				tmp[2] = nm[2];
				compute_up_dir(&(tmp[0]), render_para.updir);
				compute_up_dir(&(tmp[0]), render_para.updir);

				//do the rendering

#ifdef _SAVE_OFFRENDER_TO_IMG
				rdr.suffix = "wt_pos";
#endif
				rdr.RenderModelAt(-1, IDmask[0],
					image_stack_wt[0], temp_buffer,
					render_para.campos_pos,
					render_para.target,
					render_para.updir, sampler, cc_in_current_scene,
					OffScreenRenderer::ORTHOGRAPHIC,
					ccborder, zNear, zFar, CutEPS);

#ifdef _SAVE_OFFRENDER_TO_IMG
				rdr.suffix = "wt_neg";
#endif
				rdr.RenderModelAt(-1, IDmask[1],
					image_stack_wt[1], temp_buffer,
					render_para.campos_neg,
					render_para.target,
					render_para.updir, sampler, cc_in_current_scene,
					OffScreenRenderer::ORTHOGRAPHIC,
					ccborder, zNear, zFar, CutEPS);

#ifdef _SAVE_OFFRENDER_TO_IMG
				rdr.suffix = "wo_pos";
#endif
				rdr.RenderModelAt(i, IDmask[0],
					image_stack_wo[0], temp_buffer,
					render_para.campos_pos,
					render_para.target,
					render_para.updir,
					sampler, cc_in_current_scene,
					OffScreenRenderer::ORTHOGRAPHIC,
					ccborder, zNear, zFar, CutEPS);
#ifdef _SAVE_OFFRENDER_TO_IMG
				rdr.suffix = "wo_neg";
#endif
				rdr.RenderModelAt(i, IDmask[1],
					image_stack_wo[1], temp_buffer,
					render_para.campos_neg,
					render_para.target,
					render_para.updir,
					sampler, cc_in_current_scene,
					OffScreenRenderer::ORTHOGRAPHIC,
					ccborder, zNear, zFar, CutEPS);

				for (int batch_i = 0; batch_i < image_stack_wt.size(); batch_i++) {
					compute_the_diff_non_blank(image_stack_wo[batch_i], image_stack_wt[batch_i],
						IDmask[batch_i], i,
						1.0f, score_current_cc, total_pix_siz, img_wid);
				}

			}


			if (total_pix_siz == 0)score_current_cc = 0;
			else score_current_cc = score_current_cc / (float)total_pix_siz;
			score_list.push_back(std::make_pair(i, score_current_cc));


			//construct relation list
			for (unsigned int s = 0; s < cc_in_current_scene.size(); s++) {
				if (cc_in_current_scene[s] == false || s == i)continue;
				relation_list[s].insert(i);
			}


		}
	}
	else {
		for (auto it = update_list.begin(); it != update_list.end(); it++) {
			//ccid i
			unsigned int ccid = (*it);

			if (delete_list.find(ccid) != delete_list.end())continue;

			rdr.make_cc_NOTshow_for_compare(ccid);

			float ccborder = get_cc_border(sampler->cc_bbx[ccid]);
			ccborder *= _BORDER_TO_BBX;
			if (minimal_border>0 && ccborder < minimal_border)ccborder = minimal_border;


			float zNear = 0.1;

#ifdef _SAVE_OFFRENDER_TO_IMG
			rdr.output_counter = ccid;
			rdr.rder_couter = 0;
#endif


			std::vector<bool> cc_in_current_scene;
			cc_in_current_scene.resize(sampler->cc_sampled_triid.size(), false);
			float score_current_cc = 0.0f;
			int total_pix_siz = 0;

			for (int j = 0; j <  sampler->cc_sampled_triid[ccid].size(); j++) {
				//for each triangle in cci

				//assemble render para
				RenderStrc render_para;

				int triid = sampler->cc_sampled_triid[ccid][j].first;
				auto nm = msh->GetTriangles()[triid].normal;

				render_para.campos_pos[0] = offset_scale*nm[0] +  sampler->cc_sampled_triid[ccid][j].second[0];
				render_para.campos_pos[1] = offset_scale*nm[1] +  sampler->cc_sampled_triid[ccid][j].second[1];
				render_para.campos_pos[2] = offset_scale*nm[2] +  sampler->cc_sampled_triid[ccid][j].second[2];

				render_para.campos_neg[0] = -offset_scale*nm[0] +  sampler->cc_sampled_triid[ccid][j].second[0];
				render_para.campos_neg[1] = -offset_scale*nm[1] +  sampler->cc_sampled_triid[ccid][j].second[1];
				render_para.campos_neg[2] = -offset_scale*nm[2] +  sampler->cc_sampled_triid[ccid][j].second[2];

				render_para.target[0] =  sampler->cc_sampled_triid[ccid][j].second[0];
				render_para.target[1] =  sampler->cc_sampled_triid[ccid][j].second[1];
				render_para.target[2] =  sampler->cc_sampled_triid[ccid][j].second[2];

				std::vector<float> tmp(3);
				tmp[0] = nm[0];
				tmp[1] = nm[1];
				tmp[2] = nm[2];
				compute_up_dir(&(tmp[0]), render_para.updir);
				compute_up_dir(&(tmp[0]), render_para.updir);

				//do the rendering

#ifdef _SAVE_OFFRENDER_TO_IMG
				rdr.suffix = "wt_pos";
#endif
				rdr.RenderModelAt(-1, IDmask[0],
					image_stack_wt[0], temp_buffer,
					render_para.campos_pos,
					render_para.target,
					render_para.updir, sampler, cc_in_current_scene,
					OffScreenRenderer::ORTHOGRAPHIC,
					ccborder, zNear, zFar, CutEPS);

#ifdef _SAVE_OFFRENDER_TO_IMG
				rdr.suffix = "wt_neg";
#endif
				rdr.RenderModelAt(-1, IDmask[1],
					image_stack_wt[1], temp_buffer,
					render_para.campos_neg,
					render_para.target,
					render_para.updir, sampler, cc_in_current_scene,
					OffScreenRenderer::ORTHOGRAPHIC,
					ccborder, zNear, zFar, CutEPS);

#ifdef _SAVE_OFFRENDER_TO_IMG
				rdr.suffix = "wo_pos";
#endif
				rdr.RenderModelAt(ccid, IDmask[0],
					image_stack_wo[0], temp_buffer,
					render_para.campos_pos,
					render_para.target,
					render_para.updir,
					sampler, cc_in_current_scene,
					OffScreenRenderer::ORTHOGRAPHIC,
					ccborder, zNear, zFar, CutEPS);
#ifdef _SAVE_OFFRENDER_TO_IMG
				rdr.suffix = "wo_neg";
#endif
				rdr.RenderModelAt(ccid, IDmask[1],
					image_stack_wo[1], temp_buffer,
					render_para.campos_neg,
					render_para.target,
					render_para.updir,
					sampler, cc_in_current_scene,
					OffScreenRenderer::ORTHOGRAPHIC,
					ccborder, zNear, zFar, CutEPS);

				for (int batch_i = 0; batch_i < image_stack_wt.size(); batch_i++) {
					compute_the_diff_non_blank(image_stack_wo[batch_i], image_stack_wt[batch_i],
						IDmask[batch_i], ccid,
						1.0f, score_current_cc, total_pix_siz, img_wid);
				}

			}


			//update score
			if (total_pix_siz == 0)score_current_cc = 0;
			else score_current_cc = score_current_cc / (float)total_pix_siz;
			for (int i = 0; i < score_list.size(); i++) {
				if (score_list[i].first == ccid) {
					score_list[i].second = score_current_cc;
					break;
				}
			}

			//update relation list
			for (unsigned int s = 0; s < cc_in_current_scene.size(); s++) {
				if (cc_in_current_scene[s] == false || s == ccid)continue;
				if (delete_list.find(s) != delete_list.end()) {
					std::cout << "Weird";
#if defined(_RECORD_ERROR_LOGS)
					(*_log_off) << _file_under_processing << "ERROR preprocessing error 0\n";
					(*_log_off).flush();
#endif
				}
				relation_list[s].insert(ccid);
			}
		}

	}
}

int VisualProcesser::delete_duplicate_cc_visual_based(float offset_scale,
	float IMG_EPS, float zFar, float minimal_border, float CutEPS, std::vector<bool> &cc_keep_list) {
	OffScreenRenderer rdr(128,128);
	rdr.setMesh(msh);
	//rdr.setFrame(FRM_SIZ, FRM_SIZ);

	std::vector<std::vector<int>> IDmask;
	std::vector<int> temp_buffer;
	IDmask.resize(2, std::vector<int>(rdr.getPixelNum(), 0));
	temp_buffer.resize(rdr.getPixelNum(), 0);


	std::vector<std::vector<float>> image_stack_wt;
	image_stack_wt.resize(2, std::vector<float>(rdr.getPixelNum(), 0.0));
	std::vector<std::vector<float>> image_stack_wo;
	image_stack_wo.resize(2, std::vector<float>(rdr.getPixelNum(), 0.0));

	std::set<unsigned int> delete_list;
	std::set<unsigned int> update_list;
	std::vector<std::pair<unsigned int, float>> score_list;
	std::vector<std::set<unsigned int>> relation_list;

	bool finish_flag = false;
	int iter = 0;
	do {
		std::cout << "------------Iteration " << iter << "----------------\n";
		
		//int mborder = iter < 2 ? -1 : minimal_border;
		update_score_and_relation(rdr, offset_scale, zFar, minimal_border,  CutEPS, update_list,
			delete_list, score_list, relation_list,
			IDmask, temp_buffer, image_stack_wt, image_stack_wo);
		iter++;

		

		//sort the score list 
		std::sort(score_list.begin(), score_list.end(), [](std::pair<unsigned int, float> &a, std::pair<unsigned int, float> &b) {
			return a.second > b.second;
		});//heap?

		//gather cc that larger than eps
		std::set<int> new_exclude_ccid;

		
		std::vector<std::pair<unsigned int, float>>::iterator bkpt = score_list.end();
		for (auto it=score_list.begin(); it!=score_list.end(); it++) {

			
				//std::cout <<"debug use"<< it->second << "\n";
			
			if (it->second < IMG_EPS) {
				bkpt = it;
				break;
			}
			else {
				new_exclude_ccid.insert(it->first);
			}
		}
		//clear the relation list
		
		std::vector<int> tmp1(new_exclude_ccid.begin(), new_exclude_ccid.end());
		std::sort(tmp1.begin(), tmp1.end());
		for (int ss = 0; ss < relation_list.size(); ss++) {
			std::vector<int> v(relation_list[ss].size());
			std::vector<int> tmp2(relation_list[ss].begin(), relation_list[ss].end());
			std::sort(tmp2.begin(), tmp2.end());
			auto  it = std::set_difference(tmp2.begin(), 
				tmp2.end(), tmp1.begin(), tmp1.end(),
				v.begin());
			v.resize(it - v.begin());

			relation_list[ss].clear();
			relation_list[ss].insert(v.begin(), v.end());
		}
		//clear the score list
		if (bkpt == score_list.end()) {
			score_list.clear();
		}else if (bkpt!= score_list.begin()){
			score_list.erase(score_list.begin(), bkpt);
		}


		do {
			int siz = score_list.size();

			if (siz==0||score_list[siz - 1].second > IMG_EPS) {
				//the smallest score is already larger than the EPS
				finish_flag = true;
				break;
			}

			//parallel
			update_list.clear();
			std::vector<pair<int,unsigned int>> going_to_delete;
			int num=get_parallel_delete_list(score_list, IMG_EPS, relation_list, delete_list, going_to_delete);
			std::cout << num << " ccs will be deleted.\n";
			for (int ss = 0; ss < going_to_delete.size(); ss++) {
				int ccid = going_to_delete[ss].second;
				int score_id = going_to_delete[ss].first;
				if (ccid != score_list[score_id].first) {
					std::cout << "Error Mismatch\n";
				}
				std::cout << "CC " << ccid << " deleted(F"
					<< sampler->cc_triid[ccid][0] << "),score "
					<< score_list[score_id].second << "\n";
				delete_list.insert(ccid);
				update_list.insert(relation_list[ccid].begin(), relation_list[ccid].end());
				score_list.erase(score_list.begin() + score_id);
			}
			rdr.make_cc_NOTshow(delete_list);

		} while (update_list.size() == 0 && score_list.size() > 0);
		//relation list

		if (finish_flag)
			break;
	} while (score_list.size() > 0);

	cc_keep_list.clear();
	cc_keep_list.resize(msh->GetNumComponents(), true);
	for (auto it = delete_list.begin(); it != delete_list.end(); it++) {
		cc_keep_list[*it] = false;
	}

	return delete_list.size();
}

void VisualProcesser::save_res_mesh(const std::vector<bool> &cc_should_keep, int iter) {
	//msh->scale_back();

	if (iter == -1) {

	ofstream off("_debug_view_delete_res.obj");
	ofstream off2("_debug_view_delete_thrown.obj");
	for (int i = 0; i < msh->GetNumVertices(); i++) {
		auto vh = msh->GetVertices()[i];
		off << "v " << vh.pos[0] << " " << vh.pos[1] << " " << vh.pos[2] << " \n";
		off2 << "v " << vh.pos[0] << " " << vh.pos[1] << " " << vh.pos[2] << " \n";
	}

	for (int i = 0; i < msh->GetNumTriangles(); i++) {
		auto th = msh->GetTriangles()[i];
		if (cc_should_keep[th.component_id] == false) {
			off2 << "f " << th.vertices[0] + 1 << " " << th.vertices[1] + 1 << " " << th.vertices[2] + 1 << " \n";
		}
		else {
			off << "f " << th.vertices[0] + 1 << " " << th.vertices[1] + 1 << " " << th.vertices[2] + 1 << " \n";
		}


	}
	off2.close();
	off.close();
	}
	else {
		char buffer[128];
		sprintf(buffer, "_debug_view_delete_res_iter%d.obj", iter);
		ofstream off(buffer);
		sprintf(buffer, "_debug_view_delete_thrown_iter%d.obj", iter);
		ofstream off2(buffer);
		for (int i = 0; i < msh->GetNumVertices(); i++) {
			auto vh = msh->GetVertices()[i];
			off << "v " << vh.pos[0] << " " << vh.pos[1] << " " << vh.pos[2] << " \n";
			off2 << "v " << vh.pos[0] << " " << vh.pos[1] << " " << vh.pos[2] << " \n";
		}

		for (int i = 0; i < msh->GetNumTriangles(); i++) {
			auto th = msh->GetTriangles()[i];
			if (cc_should_keep[th.component_id] == false) {
				off2 << "f " << th.vertices[0] + 1 << " " << th.vertices[1] + 1 << " " << th.vertices[2] + 1 << " \n";
			}
			else {
				off << "f " << th.vertices[0] + 1 << " " << th.vertices[1] + 1 << " " << th.vertices[2] + 1 << " \n";
			}


		}
		off2.close();
		off.close();
	}
}


void VisualProcesser::get_unary_score_on_campos(OffScreenRenderer &rdr, float * camera, float * target, float * updir,
	std::vector<int>& posPixel, std::vector<int>& negPixel)
{
	//the triangle id in the offscreen render is from 1
	std::vector<int> IDmask;
	IDmask.resize(rdr.getPixelNum(), 0);
	rdr.RenderModelAt_scoring(camera, target, updir, 0.7f, 0.12f, 2.0f, IDmask);

	for (int i = 0; i < IDmask.size(); i++) {
		if (IDmask[i] > 0) {
			int triid = IDmask[i] - 1;

			auto nm = msh->GetTriangles()[triid].normal;

			if (abs(nm[0]) < 1e-7&&abs(nm[1]) < 1e-7&&abs(nm[2]) < 1e-7) continue;

			float dotprod = nm[0] * (target[0] - camera[0])
				+ nm[1] * (target[1] - camera[1])
				+ nm[2] * (target[2] - camera[2]);

			if (dotprod > 0) {
				negPixel[triid]++;
			}
			else {
				posPixel[triid]++;
			}
		}
	}
}

void VisualProcesser::get_unary_score_all_campos(std::vector<int> &posPix, std::vector<int> &negPix)
{
	//std::vector<int> posPix, negPix;
	posPix.clear();
	negPix.clear();
	posPix.resize(msh->GetNumTriangles(), 0);
	negPix.resize(msh->GetNumTriangles(), 0);

	msh->update_mesh_properties();

	OffScreenRenderer rdr(512,512,OffScreenRenderer::UNARY_SCORING);
	rdr.setMesh(msh);


	std::vector<vector<float>> camPoses = 
		{{0.000000f, 0.500000f, 0.000000f},
		{ 0.250000f, 0.404509f, 0.154508f },
		{ 0.250000f, 0.404509f, -0.154508f },
		{ 0.000000f, 0.425325f, -0.262866f },
		{ 0.000000f, 0.425325f, 0.262866f },
		{ 0.425325f, 0.262866f, 0.000000f },
		{ -0.250000f, 0.404509f, -0.154508f },
		{ -0.250000f, 0.404509f, 0.154508f },
		{ -0.425325f, 0.262866f, 0.000000f },
		{ 0.000000f, -0.500000f, 0.000000f },
		{ 0.250000f, -0.404509f, -0.154508f },
		{ 0.250000f, -0.404509f, 0.154508f },
		{ 0.000000f, -0.425325f, 0.262866f },
		{ 0.000000f, -0.425325f, -0.262866f },
		{ 0.425325f, -0.262866f, 0.000000f },
		{ -0.250000f, -0.404509f, 0.154508f },
		{ -0.250000f, -0.404509f, -0.154508f },
		{ -0.425325f, -0.262866f, 0.000000f },
		{ 0.500000f, 0.000000f, 0.000000f },
		{ 0.404509f, -0.154508f, -0.250000f },
		{ 0.404509f, 0.154508f, -0.250000f },
		{ 0.262866f, 0.000000f, -0.425325f } ,
		{ 0.404509f, 0.154508f, 0.250000f },
		{ 0.404509f, -0.154508f, 0.250000f },
		{ 0.262866f, 0.000000f, 0.425325f },
		{ -0.500000f, 0.000000f, 0.000000f },
		{ -0.404509f, -0.154508f, 0.250000f },
		{ -0.404509f, 0.154508f, 0.250000f },
		{ -0.262866f, 0.000000f, 0.425325f },
		{ -0.404509f, 0.154508f, -0.250000f },
		{ -0.404509f, -0.154508f, -0.250000f },
		{ -0.262866f, 0.000000f, -0.425325f },
		{ 0.000000f, 0.000000f, 0.500000f },
		{ -0.154508f, - 0.250000f, 0.404509f },//
		{ 0.154508f, -0.250000f, 0.404509f },
		{ 0.154508f, 0.250000f, 0.404509f },
		{ -0.154508f, 0.250000f, 0.404509f },
		{ 0.000000f, 0.000000f, -0.500000f },
		{ -0.154508f, 0.250000f, -0.404509f },
		{ 0.154508f, 0.250000f, -0.404509f },
		{ 0.154508f, -0.250000f, -0.404509f },
		{ -0.154508f, -0.250000f, -0.404509f } };

	for (int i = 0; i < camPoses.size(); i++) {
		float campp[3] = { camPoses[i][0] * 2.0,camPoses[i][1] * 2.0 ,camPoses[i][2] * 2.0 };
		float target[3] = { 0,0,0 };
		float updir[3];
		float camnm = sqrt(pow(campp[0], 2)
			+ pow(campp[1], 2)
			+ pow(campp[2], 2));
		float viewdir[3] = { -campp[0]/ camnm,
			-campp[1]/ camnm, -campp[2]/ camnm };
		compute_up_dir(&(viewdir[0]), updir);
		get_unary_score_on_campos(rdr, &(campp[0]), target, updir, posPix, negPix);
	}

	//std::cout << "done";
}

int VisualProcesser::get_parallel_delete_list(std::vector<std::pair<unsigned int, float>>& score_list, double DEL_EPS, 
	std::vector<std::set<unsigned int>>& relation_list,
	std::set<unsigned int>& delete_list, 
	std::vector<pair<int, unsigned int>>& need_to_del)
{
	need_to_del.clear();
	std::vector<BlackMesh::BlackMesh<double>::BBX> bbx_del_list;

	int siz = score_list.size();


	for (int i = siz - 1; i >= 0&& score_list[i].second<DEL_EPS; i--) {
		int ccid_pending = score_list[i].first;
		
		bool pass_flag = true;
		for (int j = 0; j < need_to_del.size(); j++) {
			int ccid_in = need_to_del[j].second;
			if (ccid_in == ccid_pending) {
				std::cout << "Error";
			}
			else {
				auto flag1 = if_they_have_relation(ccid_in, ccid_pending, relation_list);
				if (flag1) {
					pass_flag = false;
					break;
				}
				else {
					auto flag2 = sampler->cc_bbx[ccid_in]. if_two_bbx_intersects( sampler->cc_bbx[ccid_pending]);
					if (flag2) {
						pass_flag = false;
						break;
					}
				}
			}
		}
		
		if (pass_flag) {
			need_to_del.push_back(std::make_pair(i,ccid_pending));
			bbx_del_list.push_back(sampler->cc_bbx[ccid_pending]);
		}
	
	}

	return need_to_del.size();
}


