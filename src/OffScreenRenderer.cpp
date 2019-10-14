#include "OffScreenRenderer.h"
#include "stdafx.h"

void OffScreenRenderer::RenderModelAt(int ccid, std::vector<int> &IDmask, std::vector<float>& output_depth_img, std::vector<int> &temp_int_buf, float *pos, float *target, float * updir,
	MeshSampling *sampler, std::vector<bool> &stacked_cc, PROJECTION_TYPE type, const float border,
	const float zNear, const float zFar, const float CUTTING_EPS) {


#ifdef _USE_ZNEAR_CUT_SCENE
	//assume the target-pos dir is on normal direction
	float dist = sqrt(pow(pos[0] - target[0], 2) + pow(pos[1] - target[1], 2) + pow(pos[2] - target[2], 2));
	float newzNear = dist - CUTTING_EPS;
#else
	float newzNear = zNear;
#endif

	setCamera(pos, target, updir);
	setProjection(type, -border, border, newzNear, zFar);
	setIdentityModel();


	glBindFramebuffer(GL_FRAMEBUFFER, fbo);

	glViewport(0, 0, _width, _height);
	glClearDepth(1);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

	_shader->use();

	updateMVP();

	_shader->setMat4("MVP", _MVP);



	auto MVPInv = glm::inverse(_MVP);
	glm::vec4 p0 = MVPInv*glm::vec4(1, 1, 1, 1.0);
	glm::vec4 p1 = MVPInv*glm::vec4(1, -1, 1, 1.0);
	glm::vec4 p2 = MVPInv*glm::vec4(-1, 1, 1, 1.0);
	glm::vec4 p3 = MVPInv*glm::vec4(-1, -1, 1, 1.0);
	glm::vec4 p4 = MVPInv*glm::vec4(1, 1, -1, 1.0);
	glm::vec4 p5 = MVPInv*glm::vec4(1, -1, -1, 1.0);
	glm::vec4 p6 = MVPInv*glm::vec4(-1, 1, -1, 1.0);
	glm::vec4 p7 = MVPInv*glm::vec4(-1, -1, -1, 1.0);


	BlackMesh::BlackMesh<double>::BBX mvpbbx;
	mvpbbx.update_bbx(p0.x, p0.y, p0.z);
	mvpbbx.update_bbx(p1.x, p1.y, p1.z);
	mvpbbx.update_bbx(p2.x, p2.y, p2.z);
	mvpbbx.update_bbx(p3.x, p3.y, p3.z);
	mvpbbx.update_bbx(p4.x, p4.y, p4.z);
	mvpbbx.update_bbx(p5.x, p5.y, p5.z);
	mvpbbx.update_bbx(p6.x, p6.y, p6.z);
	mvpbbx.update_bbx(p7.x, p7.y, p7.z);


	if (ccid == -1) {
		//glBindVertexArray(_VAO);
		for (int i = 0; i < cc_visualable.size(); i++) {

			bool intflg = mvpbbx.if_two_bbx_intersects(sampler->cc_bbx[i]);
			if (cc_visualable[i]&& intflg) {
				glBindVertexArray(_VAO_ls[i]);
				glDrawArrays(GL_TRIANGLES, 0, _vtices_ls[i].size());
				glBindVertexArray(0);
			}
		}
	}
	else {

		//glBindVertexArray(_VAOC);
		make_cc_NOTshow_for_compare(ccid);
		for (int i = 0; i < cc_visualable.size(); i++) {

			bool intflg = mvpbbx.if_two_bbx_intersects(sampler->cc_bbx[i]);
			if (cc_visualable[i] && i != compare_ccidC &&(intflg)) {

				glBindVertexArray(_VAO_ls[i]);
				glDrawArrays(GL_TRIANGLES, 0, _vtices_ls[i].size());
				glBindVertexArray(0);
			}
		}
	}
	

	_shader->release();


#ifdef _SAVE_OFFRENDER_TO_IMG
	cv::Mat res_img = cv::Mat(_height, _width, CV_8UC3);
	//std::vector<std::uint8_t> data(_width*_height * 4);
	//std::vector<float> data2(_width*_height * 4);
	glReadBuffer(GL_COLOR_ATTACHMENT0);
	//glReadPixels(0, 0, _width, _height, GL_BGRA, GL_UNSIGNED_BYTE, &(data[0]));
	//glReadPixels(0, 0, _width, _height, GL_BGRA, GL_FLOAT, &(data2[0]));


	//use fast 4-byte alignment (default anyway) if possible
	glPixelStorei(GL_PACK_ALIGNMENT, (res_img.step & 3) ? 1 : 4);
	//set length of one complete row in destination data (doesn't need to equal img.cols)
	glPixelStorei(GL_PACK_ROW_LENGTH, res_img.step / res_img.elemSize());
	glReadPixels(0, 0, res_img.cols, res_img.rows, GL_BGR, GL_UNSIGNED_BYTE, res_img.data);
	cv::flip(res_img, res_img, 0);
	//glReadPixels(0, 0, width, height, GL_BGRA, GL_FLOAT, res_img.data);
	cv::imwrite("_debug_color_buffer.png", res_img);
#endif

	//glReadBuffer(GL_DEPTH_COMPONENT);//no need to call
	//std::vector<float> depth_img(_width*_height,0.0);
	//output_depth_img.resize(_width*_height, 0.0);
	glReadPixels(0, 0, _width, _height, GL_DEPTH_COMPONENT, GL_FLOAT, &(output_depth_img[0]));

	if (ccid != -1) {//deleted the cci
					 //std::vector<int> index_data(_width*_height, 0);
		glReadBuffer(GL_COLOR_ATTACHMENT1);
		glReadPixels(0, 0, _width, _height, GL_RED_INTEGER, GL_INT, &(temp_int_buf[0]));


		int total_threads = _num_threads;
		std::vector<std::vector<int>> tmp_id;
		tmp_id.resize(total_threads, std::vector<int>{});
#pragma omp parallel for
		for (int i = 0; i < _width; i++) {
			for (int j = 0; j < _width; j++) {
				int eid = i* _width + j;
				if (temp_int_buf[eid] != 0) {
					//stacked_cc[temp_int_buf[i] - 1] = true;
					int tid = omp_get_thread_num();
					tmp_id[tid].push_back(temp_int_buf[eid] - 1);
				}
			}
		}


		for (int i = 0; i < tmp_id.size(); i++) {
			for (int j = 0; j < tmp_id[i].size(); j++) {
				stacked_cc[tmp_id[i][j]] = true;
			}
		}

	}
	else {
		glReadBuffer(GL_COLOR_ATTACHMENT1);
		glReadPixels(0, 0, _width, _height, GL_RED_INTEGER, GL_INT, &(IDmask[0]));
	}

	auto p = glGetError();
#ifdef _SAVE_OFFRENDER_TO_IMG
	cv::Mat  outImg = cv::Mat(_height, _width, CV_32FC1, &(output_depth_img[0]));
	cv::flip(outImg, outImg, 0);
	char buff[128];
	sprintf(buff, "img\\_debug_depth_buffer_cc%d_%d_%s.exr", output_counter, rder_couter, suffix.c_str());
	rder_couter++;
	cv::imwrite(buff, outImg);
#endif
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void OffScreenRenderer::RenderModelAt_scoring(float *pos, float *target, float * updir,
	const float border ,const float zNear, const float zFar,
	std::vector<int> &IDmask)
{
	if (_VBO_update_flag) {//change some cc show/unshow; change camera position do not need to do this
		assemble_vtx_data();
		clearVAO();
		clearVBO();
		init_vtx_bf();
		bind_to_VAO();

		_VBO_update_flag = false;
	}
	


	float newzNear = zNear;

	setCamera(pos, target, updir);
	setProjection(ORTHOGRAPHIC, -border, border, newzNear, zFar);
	setIdentityModel();


	glBindFramebuffer(GL_FRAMEBUFFER, fbo);

	glViewport(0, 0, _width, _height);
	glClearDepth(1);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

	_shader->use();

	updateMVP();

	
	

	_shader->setMat4("MVP", _MVP);

	for (int i = 0; i < cc_visualable.size(); i++) {
		
			glBindVertexArray(_VAO_ls[i]);
			glDrawArrays(GL_TRIANGLES, 0, _vtices_ls[i].size());
			glBindVertexArray(0);
		
	}
	_shader->release();

	glReadBuffer(GL_COLOR_ATTACHMENT1);
	glReadPixels(0, 0, _width, _height, GL_RED_INTEGER, GL_INT, &(IDmask[0]));


	glBindFramebuffer(GL_FRAMEBUFFER, 0);

}
