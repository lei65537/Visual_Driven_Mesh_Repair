#pragma once
#include "OffScreenRenderer.h"
#include "BlackMesh.h"
#include "MeshSampling.h"


class VisualProcesser
{
private:
	BlackMesh::BlackMesh<double> *msh;
	MeshSampling *sampler;

	void compute_the_diff_non_blank(std::vector<float>& img1, std::vector<float>& img2,
		std::vector<int>& IDmask, int ccid,
		const  float& clear_val, float& acc_err, int& acc_msk_siz, int image_wid);

	/***********Useful Functions************/
	inline void cross_prod(float *v1, float* v2, float* output) {
		output[0] = v1[1] * v2[2] - v1[2] * v2[1];
		output[1] = v1[2] * v2[0] - v1[0] * v2[2];
		output[2] = v1[0] * v2[1] - v1[1] * v2[0];
	};


	inline void normalize(float * v) {
		float Length=sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
		if (Length < 1e-10)return;
		float invLength = (float)1 / Length;
		v[0] *= invLength;
		v[1] *= invLength;
		v[2] *= invLength;
	};

	inline void compute_up_dir(float *forward, float *up, float EPSILON = 1e-6) {
		//view dir must be normalized
		if (fabs(forward[0]) < EPSILON && fabs(forward[2]) < EPSILON)
		{
			// forward vector is pointing +Y axis
			if (forward[1] > 0) {
				up[0] = 0;
				up[1] = 0;
				up[2] = -1;
			}
			// forward vector is pointing -Y axis
			else
			{
				up[0] = 0;
				up[1] = 0;
				up[2] = 1;
			}
		}
		// in general, up vector is straight up
		else
		{
			up[0] = 0;
			up[1] = 1;
			up[2] = 0;
		}

		float left[3];
		cross_prod(up, forward, left);
		normalize(left);

		cross_prod(forward, left, up);
		normalize(up);

	};

	inline float get_cc_border(BlackMesh::BlackMesh<double>::BBX &bbx) {
		float x = bbx.xmax - bbx.xmin;
		float y = bbx.ymax - bbx.ymin;
		float z = bbx.zmax - bbx.zmin;
		float tmp = max(x, y);
		return max(tmp, z) / 2.0f;
	};

	struct RenderStrc {
		float campos_pos[3];
		float campos_neg[3];
		float target[3];
		float updir[3];
	};
	/**************************************/

public:
	VisualProcesser() {}
	~VisualProcesser() {}

	//set mesh
	inline void setMesh(BlackMesh::BlackMesh<double> * msh_) {
		msh = msh_;
	}
	inline void setSampler(MeshSampling *spl) { sampler = spl; }
	//sample  P=max{A/(diag*diag)*N, Pfix}
	
	void update_score_and_relation(OffScreenRenderer &rdr, float offset_scale, float zFar, float minimal_border,
		float CutEPS, std::set<unsigned int> &update_list, std::set<unsigned int> &delete_list,
		std::vector<std::pair<unsigned int, float>>& score_list, std::vector<std::set<unsigned int>> &relation_list,
		std::vector<std::vector<int>> &IDmask, std::vector<int> &temp_buffer,
		std::vector<std::vector<float>> &image_stack_wt, std::vector<std::vector<float>> &image_stack_wo);

	int delete_duplicate_cc_visual_based(float offset_scale, float IMG_EPS, float zFar, float minimal_border, float CutEPS, std::vector<bool> &cc_keep_list);


	void save_res_mesh(const std::vector<bool> &cc_should_keep, int iter = -1);

	


	//unary score
	void get_unary_score_on_campos(OffScreenRenderer &rdr, float *camera, float *target, float* updir,
		std::vector<int> &posPixel, std::vector<int> &negPixel);
	void get_unary_score_all_campos(std::vector<int> &posPix, std::vector<int> &negPix);

	//parallel delete
	int get_parallel_delete_list(std::vector<std::pair<unsigned int, float>> &score_list, double DEL_EPS,
		std::vector<std::set<unsigned int>> &relation_list,
		std::set<unsigned int> &delete_list,
		std::vector<pair<int, unsigned int>> &need_to_del);
	//(sampler->bbx,) score(sorted, large->small), relation_list, deleted, output, return number going to delete
private:

	inline bool if_they_have_relation(int ccidi, int ccidj, std::vector<std::set<unsigned int>>& relation_list) {

		auto it = relation_list[ccidj].find(ccidi);
		if (it == relation_list[ccidj].end()) {
			it= relation_list[ccidi].find(ccidj);
			if (it == relation_list[ccidi].end()) {
				return false;
			}
			else {
				return true;
			}
		}
		else {
			return true;
		}
				
	}
};

