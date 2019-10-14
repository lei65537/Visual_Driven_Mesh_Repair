#pragma once

#include "BlackMesh.h"

//#include <glm/glm.hpp>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

//#include <GLFW/glfw3.h>

#include "Shader.h"

//#include <opencv2/opencv.hpp>
#include <set>

#include "MeshSampling.h"

//#define _DEFAULT_FRAME_WIDTH 128
//#define _DEFAULT_FRAME_HEIGHT 128


class OffScreenRenderer
{
	//Domus Aurea
public:
	enum PROJECTION_TYPE { PERSPECTIVE = 1, ORTHOGRAPHIC };
	enum MODULE_TYPE { UNARY_SCORING = 1, DUPLICATE_DELETING};
private:
	BlackMesh::BlackMesh<double> *_msh;
	MODULE_TYPE _mdl;

	glm::mat4 _Projection;
	glm::mat4 _View;
	glm::mat4 _Model;
	glm::mat4 _MVP;

	Shader *_shader;

	std::vector<unsigned int> _VBO_ls;//binded to vertex buffer
	std::vector<unsigned int> _VAO_ls;
	std::vector<unsigned int> _VID_ls;
	bool _VBO_update_flag;
	int compare_ccidC;


	std::vector<std::vector<float>> _vtices_ls;//should not change during each iteration
	std::vector<std::vector<int>> _vccids_ls;//should not change during each iteration



	int _width;
	int _height;

	std::vector<bool> cc_visualable;

	//gl
	inline void initGlew() {
		glewExperimental = true;
		//glfwinit();
		auto flg = glewInit();
		if (flg != GLEW_OK) {
			std::cout << "GLEW INIT FAIL\n";
#if defined(_RECORD_ERROR_LOGS)
			(*_log_off) << _file_under_processing << " ERROR GLEW INIT FAIL\n";
			(*_log_off).flush();
#endif
		}
	}

	//load shader
	inline void init_shader() {
		if(_mdl== DUPLICATE_DELETING){

			_shader = new Shader(NULL, NULL);

		}
		else {

			_shader = new Shader(NULL, NULL);
		}
	}

	//init cc visualable
	inline void init_ccVisualable() {
		cc_visualable.clear();
		if (_msh->GetNumComponents() < 1)_msh->mark_component_with_coherence();
		cc_visualable.resize(_msh->GetNumComponents(), true);
	}

	//buffer
	inline void init_vtx_bf() {
		_VBO_ls.resize(cc_visualable.size(), -1);
		_VID_ls.resize(cc_visualable.size(), -1);
		for (int i = 0; i < cc_visualable.size(); i++) {
			glGenBuffers(1, &(_VBO_ls[i]));
			glBindBuffer(GL_ARRAY_BUFFER, _VBO_ls[i]);
			glBufferData(GL_ARRAY_BUFFER, sizeof(_vtices_ls[i][0]) * _vtices_ls[i].size(), &(_vtices_ls[i][0]), GL_STATIC_DRAW);

			glGenBuffers(1, &(_VID_ls[i]));
			glBindBuffer(GL_ARRAY_BUFFER, _VID_ls[i]);
			glBufferData(GL_ARRAY_BUFFER, sizeof(_vccids_ls[i][0])* _vccids_ls[i].size(), &(_vccids_ls[i][0]), GL_STATIC_DRAW);
		}
	}


	inline void bind_to_VAO() {
		_VAO_ls.resize(cc_visualable.size(), -1);
		for (int i = 0; i < cc_visualable.size(); i++) {
			glGenVertexArrays(1, &(_VAO_ls[i]));
			glBindVertexArray((_VAO_ls[i]));
			glBindBuffer(GL_ARRAY_BUFFER, (_VBO_ls[i]));

			//pos
			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
			glEnableVertexAttribArray(0);

			//ccid
			glEnableVertexAttribArray(1);
			glBindBuffer(GL_ARRAY_BUFFER, _VID_ls[i]);
			glVertexAttribIPointer(1, 1, GL_INT, 0, (void*)0);


			glBindBuffer(GL_ARRAY_BUFFER, 0);
		}
	}


	GLuint fbo;
	GLuint rboColor;
	GLuint rboIdx;
	GLuint rboDepth;
	inline void init_framebuffer() {
		// framebuffer configuration
		// create Frame Buffer Object (FBO)
		glGenFramebuffers(1, &fbo);
		glBindFramebuffer(GL_FRAMEBUFFER, fbo);

		GLenum draw_buffs[] = { GL_COLOR_ATTACHMENT0,GL_COLOR_ATTACHMENT1 };
		glDrawBuffers(2, draw_buffs);

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glGenRenderbuffers(1, &rboColor);
		glBindRenderbuffer(GL_RENDERBUFFER, rboColor);

		glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA, _width, _height);

		glBindRenderbuffer(GL_RENDERBUFFER, 0);

		////buffer for cc index

		glGenRenderbuffers(1, &rboIdx);
		glBindRenderbuffer(GL_RENDERBUFFER, rboIdx);
		glRenderbufferStorage(GL_RENDERBUFFER, GL_R32I, _width, _height);
		glBindRenderbuffer(GL_RENDERBUFFER, 0);


		// create Render Buffer Object (RBO) for depth
		//GLuint rboDepth = 0;
		glGenRenderbuffers(1, &rboDepth);
		glBindRenderbuffer(GL_RENDERBUFFER, rboDepth);
		glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, _width, _height);
		glBindRenderbuffer(GL_RENDERBUFFER, 0);
		// create Frame Buffer Object (FBO)

		glBindFramebuffer(GL_FRAMEBUFFER, fbo);
		// attach RBO to FBO
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
			GL_RENDERBUFFER, rboColor);
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1,
			GL_RENDERBUFFER, rboIdx);
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
			GL_RENDERBUFFER, rboDepth);

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}

	//data
	inline void assemble_vtx_data() {
		_vtices_ls.clear();
		_vccids_ls.clear();

		_vtices_ls.resize(cc_visualable.size(), std::vector<float>{});
		_vccids_ls.resize(cc_visualable.size(), std::vector<int>{});

		if (cc_visualable.size() == 0) {

			for (int i = 0; i < _msh->GetNumTriangles(); i++) {
				auto th = _msh->GetTriangles()[i];
				auto p1 = _msh->GetVertices()[th.vertices[0]].pos;
				auto p2 = _msh->GetVertices()[th.vertices[1]].pos;
				auto p3 = _msh->GetVertices()[th.vertices[2]].pos;

				int ccid = th.component_id;

				_vtices_ls[ccid].push_back((float)p1[0]);
				_vtices_ls[ccid].push_back((float)p1[1]);
				_vtices_ls[ccid].push_back((float)p1[2]);
				_vtices_ls[ccid].push_back((float)p2[0]);
				_vtices_ls[ccid].push_back((float)p2[1]);
				_vtices_ls[ccid].push_back((float)p2[2]);
				_vtices_ls[ccid].push_back((float)p3[0]);
				_vtices_ls[ccid].push_back((float)p3[1]);
				_vtices_ls[ccid].push_back((float)p3[2]);

				if (_mdl == DUPLICATE_DELETING) {
					_vccids_ls[ccid].push_back((int)th.component_id+1);
					_vccids_ls[ccid].push_back((int)th.component_id+1);
					_vccids_ls[ccid].push_back((int)th.component_id+1);
				}
				else {
					_vccids_ls[ccid].push_back(i+1);
					_vccids_ls[ccid].push_back(i+1);
					_vccids_ls[ccid].push_back(i+1);
				}

			}
		}
		else {
			for (int i = 0; i < _msh->GetNumTriangles(); i++) {
				auto th = _msh->GetTriangles()[i];
				auto p1 = _msh->GetVertices()[th.vertices[0]].pos;
				auto p2 = _msh->GetVertices()[th.vertices[1]].pos;
				auto p3 = _msh->GetVertices()[th.vertices[2]].pos;

				int ccid = th.component_id;

				_vtices_ls[ccid].push_back((float)p1[0]);
				_vtices_ls[ccid].push_back((float)p1[1]);
				_vtices_ls[ccid].push_back((float)p1[2]);
				_vtices_ls[ccid].push_back((float)p2[0]);
				_vtices_ls[ccid].push_back((float)p2[1]);
				_vtices_ls[ccid].push_back((float)p2[2]);
				_vtices_ls[ccid].push_back((float)p3[0]);
				_vtices_ls[ccid].push_back((float)p3[1]);
				_vtices_ls[ccid].push_back((float)p3[2]);



				if (_mdl == DUPLICATE_DELETING) {
					_vccids_ls[ccid].push_back((int)th.component_id+1);
					_vccids_ls[ccid].push_back((int)th.component_id+1);
					_vccids_ls[ccid].push_back((int)th.component_id+1);
				}
				else {
					_vccids_ls[ccid].push_back(i+1);
					_vccids_ls[ccid].push_back(i+1);
					_vccids_ls[ccid].push_back(i+1);
				}


			}
		}
	}

	

	//camera -not modifiable
	inline void setIdentityModel() {
		_Model = glm::mat4(1.0f);
	}

	inline void updateMVP() {
		_MVP = _Projection * _View * _Model;
	}

	//clear
	inline void clearVAO() {
		for (int i = 0; i < _VAO_ls.size(); i++) {
			if (_VAO_ls[i] != -1) {
				glBindVertexArray(0);
				glDeleteVertexArrays(1, &(_VAO_ls[i]));
				_VAO_ls[i] = -1;
			}
		}
		_VAO_ls.clear();
	}


	inline void clearVBO() {
		for (int i = 0; i < _VBO_ls.size(); i++) {
			if (_VBO_ls[i] != -1) {
				glBindBuffer(GL_ARRAY_BUFFER, 0);
				glDeleteBuffers(1, &(_VBO_ls[i]));
				_VBO_ls[i] = -1;
			}
		}
		_VBO_ls.clear();
	}
	
public:
	//enum PROJECTION_TYPE { PERSPECTIVE = 1, ORTHOGRAPHIC };
	//enum MODULE_TYPE {UNARY_SCORING=1, DUPLICATE_DELTEING};
	//init
	OffScreenRenderer(int width, int height, MODULE_TYPE mdl= DUPLICATE_DELETING) {
		_mdl = mdl;
		compare_ccidC = -1;

		_width = width;
		_height = height;

		_VAO_ls.clear();
		_VBO_ls.clear();
		_VID_ls.clear();

		initGlew();
		init_shader();

		init_framebuffer();
	}

	~OffScreenRenderer() {

		glDeleteFramebuffers(1, &fbo);
		glDeleteRenderbuffers(1, &rboColor);
		glDeleteRenderbuffers(1, &rboDepth);
		glDeleteRenderbuffers(1, &rboIdx);
		delete[]_shader;
	}



	//set mesh, construct vtx data/buffer, construct face map
	inline void setMesh(BlackMesh::BlackMesh<double> *msh) {
		_msh = msh;
		init_ccVisualable();
		assemble_vtx_data();
		init_vtx_bf();
		bind_to_VAO();


		_VBO_update_flag = false;
	}

	//let cc show/not show
	inline void make_cc_show(int ccid) {
		cc_visualable[ccid] = true;
	}
	inline void make_cc_NOTshow(int ccid) {
		cc_visualable[ccid] = false;
	}
	inline void make_cc_show(std::set<unsigned int> &ccidlist) {
		if (ccidlist.size() == 0)return;
		for (auto it = ccidlist.begin(); it != ccidlist.end(); it++) {
			cc_visualable[*it] = true;
		}
	}
	inline void make_cc_NOTshow(const std::set<unsigned int> &ccidlist) {
		if (ccidlist.size() == 0)return;
		for (auto it = ccidlist.begin(); it != ccidlist.end(); it++) {
			cc_visualable[*it] = false;
		}
	}
	inline void make_cc_NOTshow_for_compare(int ccid) {
		compare_ccidC = ccid;
	}
	// camera
	inline void setCamera(float *pos, float *target, float * updir) {
		_View = glm::lookAt(
			glm::vec3(pos[0], pos[1], pos[2]), // Camera is at (4,3,3), in World Space
			glm::vec3(target[0], target[1], target[2]), // and looks at the origin
			glm::vec3(updir[0], updir[1], updir[2])  // Head is up (set to 0,-1,0 to look upside-down)
		);
	}

	inline void setProjection(PROJECTION_TYPE type, const float left_bottom, const float right_top, const float zNear, const float zFar) {
		if (type == PERSPECTIVE) {
		}
		else {
			_Projection = glm::ortho(left_bottom, right_top, left_bottom, right_top, zNear, zFar);
		}
	}

#ifdef _SAVE_OFFRENDER_TO_IMG
	int output_counter;
	int rder_couter;
	string suffix;
#endif

	//get offscreen rendering
	void RenderModelAt(int ccid, std::vector<int> &IDmask, std::vector<float>& output_depth_img, std::vector<int> &temp_int_buf, float *pos, float *target, float * updir, MeshSampling *sampler,
		std::vector<bool> &stacked_cc, PROJECTION_TYPE type = ORTHOGRAPHIC, const float border = 0.7f,
		const float zNear = 0.1f, const float zFar = 2.0f, const float CUTTING_EPS = 1e-9);

	void RenderModelAt_scoring(float *pos, float *target, float * updir,
		const float border, const float zNear, const float zFar,
		std::vector<int> &IDmask);

	inline int getPixelNum() { return _width*_height; }

	inline int getWidth() { return _width; }
};