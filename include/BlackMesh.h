#pragma once
#include "stdafx.h"
#include <assert.h> 
#include <vector>
#include <iostream>
#include <Eigen\Dense>
#include <string>
#define _DEBUG_USE
#define _OUTPUT_INFO

#include "Utils.h"


namespace BlackMesh {
	template <typename Real>
	class BlackMesh
	{
	public:
		BlackMesh():mNumVertices(0), mNumEdges(0), mNumTriangles(0), mAveEdge(0), mNumComponents(0){ bbx.clear(); }
		~BlackMesh() { clear(); }

		class Vertex {
		public:
			std::vector<Real> pos;
			std::vector<int> incident_edges;
		};

		class Edge {
		public:
			//not directed
			std::pair<int,int> vertices;//first should < second
			std::vector<int> incident_faces;
			bool is_nonmanifold_edge=-1;

		};

		class Triangle {
		public:
			std::vector<int> vertices;
			std::vector<int> incident_edges;
			std::vector<double> normal;
			std::vector<double> face_center;//true face center
			int component_id;
		};

		class BBX {
		private:
			inline bool if_two_bbx_intersects_axis(Real xmin1, Real xmax1, Real xmin2, Real xmax2) {
				if (xmax1 < xmin2 || xmax2 < xmin1) {
					return false;
				}
				else {
					return true;
				}
			}
		public:
			Real xmin;
			Real xmax;
			Real ymin;
			Real ymax;
			Real zmin;
			Real zmax;
			Real diag_dist;
			BBX() {
				xmin = 10000;
				xmax = -10000;
				ymin = 10000;
				ymax = -10000;
				zmin = 10000;
				zmax = -10000;
				diag_dist = 0;
				
			}
			inline void clear() {
				xmin = 10000;
				xmax = -10000;
				ymin = 10000;
				ymax = -10000;
				zmin = 10000;
				zmax = -10000;
				diag_dist = 0;

			}
			inline Real update_diag_dist() {
				diag_dist = sqrt(pow(xmax - xmin, 2) + pow(ymax - ymin, 2) + pow(zmax - zmin, 2));
				return diag_dist;
			}

			inline void update_bbx(Real* pos) {
				if (xmin > pos[0]) {
					xmin = pos[0];
				}
				if (ymin > pos[1]) {
					ymin = pos[1];
				}
				if (zmin > pos[2]) {
					zmin = pos[2];
				}
				if (xmax < pos[0]) {
					xmax = pos[0];
				}
				if (ymax < pos[1]) {
					ymax = pos[1];
				}
				if (zmax < pos[2]) {
					zmax = pos[2];
				}
			}

			inline void update_bbx(Real pos0, Real pos1, Real pos2) {
				if (xmin > pos0) {
					xmin = pos0;
				}
				if (ymin > pos1) {
					ymin = pos1;
				}
				if (zmin > pos2) {
					zmin = pos2;
				}
				if (xmax < pos0) {
					xmax = pos0;
				}
				if (ymax < pos1) {
					ymax = pos1;
				}
				if (zmax < pos2) {
					zmax = pos2;
				}
			}

			inline bool if_two_bbx_intersects(BBX &bbx) {
				if (!if_two_bbx_intersects_axis(xmin, xmax, bbx.xmin, bbx.xmax)) {
					return false;
				}
				else if (!if_two_bbx_intersects_axis(ymin, ymax, bbx.ymin, bbx.ymax)) {
					return false;
				}
				else if (!if_two_bbx_intersects_axis(zmin, zmax, bbx.zmin, bbx.zmax)) {
					return false;
				}
				else {
					return true;
				}
			}
		};

	protected:
		int mNumVertices, mNumEdges, mNumTriangles;
		int mNumComponents;
		std::vector<Vertex> mVertices;
		std::vector<Edge> mEdges;
		std::vector<Triangle> mTriangles;

		Real mAveEdge;
		int max_edge_degree;

		int if_edge_already_exist(int i, int j)//;//return id or -1
		{
			assert(i != j);
			int query_id_min, query_id_max;
			if (i > j) { query_id_max = i; query_id_min = j; }
			else { query_id_max = j; query_id_min = i; }


			for (int k = 0; k < mVertices[i].incident_edges.size(); k++) {
				int edgeid = mVertices[i].incident_edges[k];
				Edge tmp = mEdges[edgeid];
				if (tmp.vertices.first == query_id_min&&tmp.vertices.second == query_id_max) {
					return edgeid;
				}

			}

			return -1;
		}

		Real scale_factor;
		std::vector<Real> translation_factor;

	public:
		//bounding box
		BBX bbx;

		inline void clear() {
			mNumVertices = mNumEdges = mNumTriangles= mNumComponents = 0;
			mVertices.clear();
			mEdges.clear();
			mTriangles.clear();

			mAveEdge = 0;
			max_edge_degree = 0;

			stacked_V.resize(0, 0);
			stacked_F.resize(0, 0);
			stacked_face_center.resize(0, 0);
			face_color_.resize(0, 0);
			bbx.clear();
			scale_factor = 0;
			translation_factor.clear();
		}

		inline void clear_but_keep_vtx() {
			mNumEdges = mNumTriangles= mNumComponents = 0;
			mEdges.clear();
			mTriangles.clear();

			max_edge_degree = 0; 
			mAveEdge = 0;

			stacked_F.resize(0, 0);
			stacked_face_center.resize(0, 0);
			face_color_.resize(0, 0);

			for (int i = 0; i < mVertices.size(); i++) {
				mVertices[i].incident_edges.clear();
			}
		}

		inline void clear_component_id() {
			mNumComponents = 0;
			for (int i = 0; i < mTriangles.size(); i++) {
				mTriangles[i].component_id=-1;
			}
		}

		inline double GetAverageEdge() const { return Utils::to_double(mAveEdge); }
		 inline int GetNumVertices() const { return mNumVertices; }
		 inline int GetNumEdges() const { return mNumEdges; }
		 inline int GetNumTriangles() const { return mNumTriangles; }
		 inline int GetNumComponents() const { return mNumComponents; }
		 inline int GetMaxEdgeDegree() const { return max_edge_degree; }
		 inline const Vertex* GetVertices() const { return &(mVertices[0]); }
		 inline Vertex* GetVerticesEditable() { return &(mVertices[0]); }
		 inline const std::vector<Vertex> GetVerticesCopy() const { return mVertices; }
		 inline const Edge* GetEdges() const { return &(mEdges[0]); }
		 inline Edge* GetEdgesEditable()  { return &(mEdges[0]); }
		 inline const Triangle* GetTriangles() const { return &(mTriangles[0]); }
		 inline Triangle* GetTrianglesEditable() { return &(mTriangles[0]); }
		 inline const std::vector<Triangle> GetTrianglesCopy() const { return mTriangles; }
		 inline void clearNumComponent() { mNumComponents = 0; }

		void insert_vtx(std::vector<Real>& pos)//;//just insert vtx
		{
			BlackMesh::Vertex vtx;
			vtx.pos.reserve(3);
			vtx.pos = pos;
			mVertices.push_back(vtx);
			BlackMesh::mNumVertices++;

			if (bbx.xmin > pos[0]) {
				bbx.xmin = pos[0];
			}
			if (bbx.ymin > pos[1]) {
				bbx.ymin = pos[1];
			}
			if (bbx.zmin > pos[2]) {
				bbx.zmin = pos[2];
			}
			if (bbx.xmax < pos[0]) {
				bbx.xmax = pos[0];
			}
			if (bbx.ymax < pos[1]) {
				bbx.ymax = pos[1];
			}
			if (bbx.zmax < pos[2]) {
				bbx.zmax = pos[2];
			}
		}

		void insert_face(std::vector<int>& indices)//;//insert face and insert edge, also update adjacent structure.
		{
			BlackMesh::Triangle face;
			face.vertices.reserve(3);
			face.vertices = indices;
			face.component_id = -1;

			mTriangles.push_back(face);
			mNumTriangles++;

			//insert edge
			std::vector<std::pair<int, int>> edges;
			int siz_edges = indices.size();
			for (int i = 0; i < siz_edges; i++) {
				int id0 = indices[i];
				int id1 = indices[(i + 1) % siz_edges];

				//if edge doese not show in vetex->incident vertex, create one, add in vtx adjacent list
				// otherwise return the id of edge
				int edge_exist = if_edge_already_exist(id0, id1);
				if (edge_exist == -1) {
					//edge not exist
					Edge tmp;
					int query_id_min, query_id_max;
					if (id0 > id1) { query_id_max = id0; query_id_min = id1; }
					else { query_id_max = id1; query_id_min = id0; }

					tmp.vertices.first = query_id_min;
					tmp.vertices.second = query_id_max;

					mEdges.push_back(tmp);
					mNumEdges++;

					int new_edge_id = mEdges.size();
					new_edge_id--;
					edge_exist = new_edge_id;

					mVertices[id0].incident_edges.push_back(new_edge_id);
					mVertices[id1].incident_edges.push_back(new_edge_id);
				}

				mEdges[edge_exist].incident_faces.push_back(mNumTriangles - 1);
				mTriangles[mNumTriangles - 1].incident_edges.push_back(edge_exist);

			}


		}

		void copy_vertice_from(BlackMesh *inputmesh)//;
		{
			bbx.clear();

			for (int i = 0; i < inputmesh->GetNumVertices(); i++) {
				auto pos = inputmesh->GetVertices()[i].pos;
				insert_vtx(pos);
			}
		}

		int update_edge_manifold_flag()//;
		{
			auto if_two_face_coherent = [](BlackMesh::BlackMesh<Real>::Triangle &tr1,
				BlackMesh::BlackMesh<Real>::Triangle &tr2) {
				std::vector<int> vtx1 = tr1.vertices;
				std::vector<int> vtx2 = tr2.vertices;
				for (int i = 0; i < vtx1.size(); i++) {
					for (int j = 0; j < vtx2.size(); j++) {
						if (vtx1[i] == vtx2[j]) {
							if (vtx1[(i + 1) % 3] == vtx2[(j + 2) % 3] || vtx1[(i + 2) % 3] == vtx2[(j + 1) % 3]) return true;
						}
					}
				}
				return false;
			};

			int ct = 0;
			max_edge_degree = 0;

			for (int i = 0; i < mEdges.size(); i++) {
				if (mEdges[i].incident_faces.size() > max_edge_degree)max_edge_degree = mEdges[i].incident_faces.size();

				if (mEdges[i].incident_faces.size() > 2) {
					mEdges[i].is_nonmanifold_edge = 1;
					ct++;
				}
				else {
					if (mEdges[i].incident_faces.size() == 2) {
						int id1 = mEdges[i].incident_faces[0];
						int id2 = mEdges[i].incident_faces[1];
						bool flg = if_two_face_coherent(mTriangles[id1], mTriangles[id2]);
						if (!flg) {
							mEdges[i].is_nonmanifold_edge = 1;
							ct++;
							continue;
						}
					}
					mEdges[i].is_nonmanifold_edge = 0;
				}

			}
			return ct;
		}

		void compute_normal()//;
		{
			for (int i = 0; i < mTriangles.size(); i++) {
				int id0 = mTriangles[i].vertices[0];
				int id1 = mTriangles[i].vertices[1];
				int id2 = mTriangles[i].vertices[2];

				std::vector<double> p0, p1, p2;
				Utils::to_double<Real>(mVertices[id0].pos, p0);
				Utils::to_double<Real>(mVertices[id1].pos, p1);
				Utils::to_double<Real>(mVertices[id2].pos, p2);


				Eigen::Vector3d vp0(p0.data());
				Eigen::Vector3d vp1(p1.data());
				Eigen::Vector3d vp2(p2.data());

				Eigen::Vector3d vp01, vp02, nm;
				vp01 = vp1 - vp0; vp02 = vp2 - vp0;
				nm = vp01.cross(vp02);
				nm.normalize();

				mTriangles[i].normal = std::vector<double>{ nm[0],nm[1],nm[2] };
			}
		}

		void compute_face_center()//;
		{
			for (int i = 0; i < mTriangles.size(); i++) {
				mTriangles[i].face_center = { 0, 0, 0 };
				for (int j = 0; j < mTriangles[i].vertices.size(); j++) {
					std::vector<double> pt;
					Utils::to_double<Real>(mVertices[mTriangles[i].vertices[j]].pos, pt);
					mTriangles[i].face_center[0] += pt[0];
					mTriangles[i].face_center[1] += pt[1];
					mTriangles[i].face_center[2] += pt[2];
				}
				mTriangles[i].face_center[0] = mTriangles[i].face_center[0] / (float)mTriangles[i].vertices.size();
				mTriangles[i].face_center[1] = mTriangles[i].face_center[1] / (float)mTriangles[i].vertices.size();
				mTriangles[i].face_center[2] = mTriangles[i].face_center[2] / (float)mTriangles[i].vertices.size();
			}
		}

		void compute_avg_edge()//;

		{
			double lenth = 0;
			for (int i = 0; i<mEdges.size(); i++) {

				std::vector<double> pt1, pt2;
				Utils::to_double<Real>(mVertices[mEdges[i].vertices.first].pos, pt1);
				Utils::to_double<Real>(mVertices[mEdges[i].vertices.second].pos, pt2);

				Eigen::Vector3d vp1(pt1.data());
				Eigen::Vector3d vp2(pt2.data());

				double tmp = (vp1 - vp2).norm();
				lenth += tmp;
			}

			mAveEdge = lenth / (double)mEdges.size();

		}

		void update_mesh_properties()//;

		{
			assert(mNumEdges = mEdges.size());
			assert(mNumVertices = mVertices.size());
			assert(mNumTriangles = mTriangles.size());
			compute_normal();
			compute_face_center();
			compute_avg_edge();
			assemble_stacked_V_and_F();
		}

		//a stacked version of V and F, first dimension is 3!
		Eigen::MatrixXf stacked_V;
		Eigen::MatrixXi stacked_F;
		Eigen::MatrixXf stacked_face_center;
		void assemble_stacked_V_and_F()//;

		{
			stacked_V.resize(3, mNumVertices);
			stacked_F.resize(3, mNumTriangles);
			stacked_face_center.resize(3, mNumTriangles);
			for (int i = 0; i < mNumVertices; i++) {
				stacked_V.col(i) = Eigen::Vector3f(Utils::to_double(mVertices[i].pos[0]),
					Utils::to_double(mVertices[i].pos[1]), Utils::to_double(mVertices[i].pos[2]));
			}
			for (int i = 0; i < mNumTriangles; i++) {
				stacked_F.col(i) = Eigen::Vector3i(mTriangles[i].vertices.data());
				std::vector<float> tmp = { (float)mTriangles[i].face_center[0],
					(float)mTriangles[i].face_center[1],(float)mTriangles[i].face_center[2] };
				stacked_face_center.col(i) = Eigen::Vector3f(tmp.data());
			}
		}

		MatrixXu8 face_color_;
		bool load_face_color(const std::string& _filename, std::string* _msg)//;
		{
			FILE * fp = fopen(_filename.c_str(), "r");
			if (fp == nullptr) {
				if (_msg) _msg->assign("Open color file failed!");
				return false;
			}

			float b;
			std::vector<float> vecbuffer;
			while (fscanf(fp, "%f", &b) != EOF) {
				vecbuffer.push_back(b);
			}
			fclose(fp);

			if (vecbuffer.size() == mNumTriangles) {
				face_color_.resize(3, mNumTriangles);
				LibGP::color_map(face_color_, vecbuffer, LibGP::JET);
				if (_msg) _msg->clear();
				return true;
			}
			else if (vecbuffer.size() == 3 * mNumTriangles) {
				float vmin = *std::min_element(vecbuffer.begin(), vecbuffer.end());
				float vmax = *std::max_element(vecbuffer.begin(), vecbuffer.end()) + EPS;

				face_color_.resize(3, mNumTriangles);
				for (int i = 0; i < mNumTriangles; i++) {
					face_color_.col(i) <<
						UINT8(255.0 * (vecbuffer[3 * i + 0] - vmin) / (vmax - vmin)),
						UINT8(255.0 * (vecbuffer[3 * i + 1] - vmin) / (vmax - vmin)),
						UINT8(255.0 * (vecbuffer[3 * i + 2] - vmin) / (vmax - vmin));
				}

				if (_msg) _msg->clear();
				return true;
			}
			else {
				if (_msg) _msg->assign("Face color number error!");
				return false;
			}
		}

		//mark component
		
		int mark_component()//;

		{
			mNumComponents = 0;
			clear_component_id();

			for (int i = 0; i < mNumTriangles; i++) {
				if (mTriangles[i].component_id != -1) continue;

				std::queue<int> que;
				que.push(i);

				while (que.size() != 0) {
					int id = que.front();
					que.pop();

					//mark id as visited
					mTriangles[id].component_id = mNumComponents;

					//push all incident faces unvistied to the queue
					for (int j = 0; j < mTriangles[id].incident_edges.size(); j++) {
						std::vector<int> icd_faces = mEdges[mTriangles[id].incident_edges[j]].incident_faces;
						for (int s = 0; s < icd_faces.size(); s++) {
							int potential_face = icd_faces[s];
							if (mTriangles[potential_face].component_id == -1) que.push(potential_face);
						}
					}
				}

				mNumComponents++;
			}

			return mNumComponents;
		}
		int mark_component_with_coherence()//;

		{
			mNumComponents = 0;
			update_edge_manifold_flag();
			clear_component_id();

			for (int i = 0; i < mNumTriangles; i++) {
				if (mTriangles[i].component_id != -1) continue;

				std::queue<int> que;
				que.push(i);
				//mark id as visited
				mTriangles[i].component_id = mNumComponents;

				while (que.size() != 0) {
					int id = que.front();
					que.pop();

					//push all incident faces unvistied to the queue
					for (int j = 0; j < mTriangles[id].incident_edges.size(); j++) {
						std::vector<int> icd_faces = mEdges[mTriangles[id].incident_edges[j]].incident_faces;
						for (int s = 0; s < icd_faces.size(); s++) {
							int potential_face = icd_faces[s];
							if (mTriangles[potential_face].component_id == -1 && mEdges[mTriangles[id].incident_edges[j]].is_nonmanifold_edge == 0)
							{
								que.push(potential_face);
								mTriangles[potential_face].component_id = mNumComponents;
							}
						}
					}
				}

				mNumComponents++;
			}

			return mNumComponents;
		}

		//generate face color according to component id
		void generate_face_color()//;

		{
			
			face_color_.resize(3, mNumTriangles);
			for (int i = 0; i < mNumTriangles; i++) {
				int idd = mTriangles[i].component_id % 17;
				face_color_(0, i) = (UINT8)(color_table[idd][0] * 255);
				face_color_(1, i) = (UINT8)(color_table[idd][1] * 255);
				face_color_(2, i) = (UINT8)(color_table[idd][2] * 255);
			}
		}

		//split one edge at midpoint
		void split_edge(int edgeid)//;

		{
			int vtx1 = mEdges[edgeid].vertices.first;
			int vtx2 = mEdges[edgeid].vertices.second;

			std::vector<BlackMesh::BlackMesh::Triangle> tri_copy = GetTrianglesCopy();

			clear_but_keep_vtx();

			insert_vtx(std::vector<Real>{0.5*(mVertices[vtx1].pos[0] + mVertices[vtx2].pos[0]),
				0.5*(mVertices[vtx1].pos[1] + mVertices[vtx2].pos[1]),
				0.5*(mVertices[vtx1].pos[2] + mVertices[vtx2].pos[2]) });

			int newvid = mNumVertices - 1;

			for (int i = 0; i < tri_copy.size(); i++) {



				int hitflag = 0;
				for (int j = 0; j < tri_copy[i].vertices.size(); j++) {
					if (tri_copy[i].vertices[j] == vtx1 || tri_copy[i].vertices[j] == vtx2) {
						hitflag++;
					}
				}

				if (hitflag != 2) {
					insert_face(tri_copy[i].vertices);
				}
				else {

					std::vector<int> temp;
					for (int j = 0; j < tri_copy[i].vertices.size(); j++) {

						if (tri_copy[i].vertices[j] == vtx1) {
							temp.push_back(newvid);
						}
						else {
							temp.push_back(tri_copy[i].vertices[j]);
						}

					}

					insert_face(temp);

					temp.clear();
					for (int j = 0; j < tri_copy[i].vertices.size(); j++) {

						if (tri_copy[i].vertices[j] == vtx2) {
							temp.push_back(newvid);
						}
						else {
							temp.push_back(tri_copy[i].vertices[j]);
						}

					}

					insert_face(temp);
				}
			}

			update_mesh_properties();

		}
		void split_edge(std::vector<int> &edgeidlist)//;
		{

			if (edgeidlist.size() == 0) return;

			std::map<std::pair<int, int>, int> vticesid2mid;

			std::vector<BlackMesh::BlackMesh<Real>::Triangle> tri_copy = GetTrianglesCopy();

			auto edgeCopy = mEdges;

			clear_but_keep_vtx();

			for (int i = 0; i < edgeidlist.size(); i++) {
				int vtx1 = edgeCopy[edgeidlist[i]].vertices.first;
				int vtx2 = edgeCopy[edgeidlist[i]].vertices.second;

				insert_vtx(std::vector<Real>{0.5*(mVertices[vtx1].pos[0] + mVertices[vtx2].pos[0]),
					0.5*(mVertices[vtx1].pos[1] + mVertices[vtx2].pos[1]),
					0.5*(mVertices[vtx1].pos[2] + mVertices[vtx2].pos[2]) });

				int newvid = mNumVertices - 1;

				//make sure i< j in pair<i,j>
				if (vtx1<vtx2)
					vticesid2mid.insert(std::make_pair(std::make_pair(vtx1, vtx2), newvid));
				else {
					vticesid2mid.insert(std::make_pair(std::make_pair(vtx2, vtx1), newvid));
				}
			}


			for (int i = 0; i < tri_copy.size(); i++) {
				std::vector<int> hitpos;
				for (int j = 0; j < tri_copy[i].vertices.size(); j++) {
					int vtx1 = tri_copy[i].vertices[j];
					int vtx2 = tri_copy[i].vertices[(j + 1) % 3];

					if (vtx1 < vtx2) {
						auto it = vticesid2mid.find(std::make_pair(vtx1, vtx2));
						if (it != vticesid2mid.end()) {
							hitpos.push_back(j);
						}
					}
					else {
						auto it = vticesid2mid.find(std::make_pair(vtx2, vtx1));
						if (it != vticesid2mid.end()) {
							hitpos.push_back(j);
						}
					}
				}

				if (hitpos.size() == 0) {
					insert_face(tri_copy[i].vertices);
				}
				else if (hitpos.size() == 1) {
					int vtx1 = tri_copy[i].vertices[hitpos[0]];
					int vtx2 = tri_copy[i].vertices[(hitpos[0] + 1) % 3];

					int minvtx = vtx1 < vtx2 ? vtx1 : vtx2;
					int maxvtx = vtx1 > vtx2 ? vtx1 : vtx2;
					auto it = vticesid2mid.find(std::make_pair(minvtx, maxvtx));

					insert_face(std::vector<int>{it->second,
						tri_copy[i].vertices[(hitpos[0] + 1) % 3],
						tri_copy[i].vertices[(hitpos[0] + 2) % 3]});

					insert_face(std::vector<int>{tri_copy[i].vertices[hitpos[0]],
						it->second,
						tri_copy[i].vertices[(hitpos[0] + 2) % 3]});
				}
				else if (hitpos.size() == 2) {
					int vtx11 = tri_copy[i].vertices[hitpos[0]];
					int vtx12 = tri_copy[i].vertices[(hitpos[0] + 1) % 3];
					int minvtx1 = vtx11 < vtx12 ? vtx11 : vtx12;
					int maxvtx1 = vtx11 > vtx12 ? vtx11 : vtx12;

					int vtx21 = tri_copy[i].vertices[hitpos[1]];
					int vtx22 = tri_copy[i].vertices[(hitpos[1] + 1) % 3];
					int minvtx2 = vtx21 < vtx22 ? vtx21 : vtx22;
					int maxvtx2 = vtx21 > vtx22 ? vtx21 : vtx22;
					if (vtx21 == vtx12) {
						auto it1 = vticesid2mid.find(std::make_pair(minvtx1, maxvtx1));
						auto it2 = vticesid2mid.find(std::make_pair(minvtx2, maxvtx2));

						insert_face(std::vector<int>{vtx21,
							it2->second,
							it1->second});

						insert_face(std::vector<int>{vtx22,
							it1->second,
							it2->second});

						insert_face(std::vector<int>{vtx22,
							vtx11,
							it1->second});
					}
					else if (vtx22 = vtx11) {
						auto it1 = vticesid2mid.find(std::make_pair(minvtx1, maxvtx1));
						auto it2 = vticesid2mid.find(std::make_pair(minvtx2, maxvtx2));

						insert_face(std::vector<int>{vtx11,
							it1->second,
							it2->second});

						insert_face(std::vector<int>{vtx12,
							it2->second,
							it1->second});

						insert_face(std::vector<int>{vtx12,
							vtx21,
							it2->second});
					}
				}
				else if (hitpos.size() == 3) {
					auto it0 = vticesid2mid.find(std::make_pair(min(tri_copy[i].vertices[0], tri_copy[i].vertices[1]),
						max(tri_copy[i].vertices[0], tri_copy[i].vertices[1])));
					auto it1 = vticesid2mid.find(std::make_pair(min(tri_copy[i].vertices[1], tri_copy[i].vertices[2]),
						max(tri_copy[i].vertices[1], tri_copy[i].vertices[2])));
					auto it2 = vticesid2mid.find(std::make_pair(min(tri_copy[i].vertices[0], tri_copy[i].vertices[2]),
						max(tri_copy[i].vertices[0], tri_copy[i].vertices[2])));

					insert_face(std::vector<int>{tri_copy[i].vertices[0],
						it0->second,
						it2->second});

					insert_face(std::vector<int>{tri_copy[i].vertices[1],
						it1->second,
						it0->second});

					insert_face(std::vector<int>{tri_copy[i].vertices[2],
						it2->second,
						it1->second});

					insert_face(std::vector<int>{it0->second,
						it1->second,
						it2->second});
				}
				else {
					std::cout << "This should not happen\n";
#if defined(_RECORD_ERROR_LOGS)
					(*_log_off) << _file_under_processing << " ERROR 2\n";
					(*_log_off).flush();
#endif
				}
			}

			update_mesh_properties();
		}

		//scale to unit
		void scale_to_unit()//;

		{
			Real xmid = 0.5*(bbx.xmax + bbx.xmin);
			Real ymid = 0.5*(bbx.ymax + bbx.ymin);
			Real zmid = 0.5*(bbx.zmax + bbx.zmin);

			//old mesh + trans -> new mesh
			translation_factor.clear();
			translation_factor.reserve(3);
			translation_factor.push_back(-xmid);
			translation_factor.push_back(-ymid);
			translation_factor.push_back(-zmid);

			bbx.xmax += translation_factor[0];
			bbx.ymax += translation_factor[1];
			bbx.zmax += translation_factor[2];

			bbx.xmin += translation_factor[0];
			bbx.ymin += translation_factor[1];
			bbx.zmin += translation_factor[2];

			////translated mesh *scale =new mesh
			Real tmp = max(bbx.xmax, bbx.ymax);
			tmp = max(tmp, bbx.zmax);
			scale_factor = 0.5 / tmp;

			bbx.xmax *= scale_factor;
			bbx.ymax *= scale_factor;
			bbx.zmax *= scale_factor;

			bbx.xmin *= scale_factor;
			bbx.ymin *= scale_factor;
			bbx.zmin *= scale_factor;

			for (int i = 0; i < mVertices.size(); i++) {
				mVertices[i].pos[0] = (mVertices[i].pos[0] + translation_factor[0])*scale_factor;
				mVertices[i].pos[1] = (mVertices[i].pos[1] + translation_factor[1])*scale_factor;
				mVertices[i].pos[2] = (mVertices[i].pos[2] + translation_factor[2])*scale_factor;
			}

			update_mesh_properties();
		}
		void scale_back()//; 
		{
			for (int i = 0; i < mVertices.size(); i++) {
				mVertices[i].pos[0] = (mVertices[i].pos[0] / scale_factor) - translation_factor[0];
				mVertices[i].pos[1] = (mVertices[i].pos[1] / scale_factor) - translation_factor[1];
				mVertices[i].pos[2] = (mVertices[i].pos[2] / scale_factor) - translation_factor[2];
			}

		bbx.xmax = (bbx.xmax / scale_factor) - translation_factor[0];
		bbx.ymax = (bbx.ymax / scale_factor) - translation_factor[0];
		bbx.zmax = (bbx.zmax / scale_factor) - translation_factor[0];

		bbx.xmin = (bbx.xmax / scale_factor) - translation_factor[0];
		bbx.ymin = (bbx.ymax / scale_factor) - translation_factor[0];
		bbx.zmin = (bbx.zmax / scale_factor) - translation_factor[0];

		update_mesh_properties();
		}

		//for score computing
		double compute_fake_height(int triid, int edgeid, double *height, int depth, int maxdepth)//;
		{
			double nmres[3];
			double succ = compute_fake_normal(triid, edgeid, nmres, depth, maxdepth);

			if (succ<0) return -1.0;

			int vidf = -1;
			int vide = -1;
			for (int i = 0; i < 3; i++) {
				if ((mTriangles[triid].vertices[i] == mEdges[edgeid].vertices.first
					&&mTriangles[triid].vertices[(i + 1) % 3] == mEdges[edgeid].vertices.second) ||
					(mTriangles[triid].vertices[i] == mEdges[edgeid].vertices.second
						&&mTriangles[triid].vertices[(i + 1) % 3] == mEdges[edgeid].vertices.first)) {
					vidf = mTriangles[triid].vertices[i];
					vide = mTriangles[triid].vertices[(i + 1) % 3];
				}
			}
			Eigen::Vector3d nmeig;
			nmeig << nmres[0], nmres[1], nmres[2];

			Eigen::Vector3d egeig;
			egeig << Utils::to_double(mVertices[vide].pos[0] - mVertices[vidf].pos[0]),
				Utils::to_double(mVertices[vide].pos[1] - mVertices[vidf].pos[1]),
				Utils::to_double(mVertices[vide].pos[2] - mVertices[vidf].pos[2]);

			egeig.normalize();

			Eigen::Vector3d hei = nmeig.cross(egeig);
			height[0] = hei[0];
			height[1] = hei[1];
			height[2] = hei[2];

			return succ;
		}

		double compute_fake_normal(int triid, int edgeid, double * nmres, int depth, int maxdepth)
			{
				auto getTriangleArea = [&](const int &id) {
			#define x 0
			#define y 1
			#define z 2
					auto tri = mTriangles[id];
			
					std::vector<double> p0, p1, p2;
					Utils::to_double<Real>(mVertices[tri.vertices[0]].pos, p0);
					Utils::to_double<Real>(mVertices[tri.vertices[1]].pos, p1);
					Utils::to_double<Real>(mVertices[tri.vertices[2]].pos, p2);
			
					return 0.5*abs(((p1[y] - p0[y])*(p2[z] - p0[z]) + (p1[z] - p0[z])*(p2[x] - p0[x]) + (p1[x] - p0[x])*(p2[y] - p0[y])) -
						((p2[y] - p0[y])*(p1[z] - p0[z]) + (p2[z] - p0[z])*(p1[x] - p0[x]) + (p2[x] - p0[x])*(p1[y] - p0[y])));
			#undef x
			#undef y
			#undef z
				};
			
				if (depth > maxdepth)return -1.0;
			
				auto th = mTriangles[triid];
			
				auto nm = th.normal;
			
				std::vector<double> ap;
				auto eh = mEdges[edgeid];
				int vhid1 = eh.vertices.first;
				int vhid2 = eh.vertices.second;
			
				if (abs(nm[0]) < 1e-7&&abs(nm[1]) < 1e-7&&abs(nm[2]) < 1e-7) {
					auto ehid1 = -1;
					auto ehid2 = -1;
			
					int inideh = -1;
					for (int i = 0; i < 3; i++) {
						if ((th.vertices[i] == vhid1&& th.vertices[(i + 1) % 3] == vhid2) ||
							(th.vertices[i] == vhid2&& th.vertices[(i + 1) % 3] == vhid1)) {
							inideh = i;
							break;
						}
					}
			
			#ifdef _OUTPUT_DEBUG_LOGS
					if (inideh == -1) {
						std::cout << "Error";
			#if defined(_RECORD_ERROR_LOGS)
						(*_log_off) << _file_under_processing << " ERROR 3\n";
						(*_log_off).flush();
			#endif
					}
			#endif
					int anothervtx = th.vertices[(inideh + 2) % 3];
			
					auto anothervh = mVertices[anothervtx];
					for (int i = 0; i < anothervh.incident_edges.size(); i++) {
						auto ehidtmp = anothervh.incident_edges[i];
						if (mEdges[ehidtmp].vertices.first == vhid1 || mEdges[ehidtmp].vertices.second == vhid1) {
							ehid1 = ehidtmp;
						}
						if (mEdges[ehidtmp].vertices.first == vhid2 || mEdges[ehidtmp].vertices.second == vhid2) {
							ehid2 = ehidtmp;
						}
						if (ehid1 != -1 && ehid2 != -1)break;
					}
			
			#ifdef _OUTPUT_DEBUG_LOGS
					if (ehid1 == -1 || ehid2 == -1) {
						std::cout << "Error";
			#if defined(_RECORD_ERROR_LOGS)
						(*_log_off) << _file_under_processing << " ERROR 4\n";
						(*_log_off).flush();
			#endif
					}
			#endif
			
					//if the two edges are not non manifold, fetch its vtx
					double flag = -1;
					if (mEdges[ehid1].incident_faces.size() == 2) {
						int nexttriid = -1;
						nexttriid = mEdges[ehid1].incident_faces[0] == triid ?
							mEdges[ehid1].incident_faces[1] : mEdges[ehid1].incident_faces[0];
						flag = compute_fake_height(nexttriid, ehid1, nmres, depth + 1, maxdepth);
					}
					if (flag>0) return flag;
			
					if (mEdges[ehid2].incident_faces.size() == 2) {
						int nexttriid = -1;
						nexttriid = mEdges[ehid2].incident_faces[0] == triid ?
							mEdges[ehid2].incident_faces[1] : mEdges[ehid2].incident_faces[0];
						flag = compute_fake_height(nexttriid, ehid2, nmres, depth + 1, maxdepth);
					}
					if (flag>0) return flag;
				}
				else {
			
					nmres[0] = nm[0];
					nmres[1] = nm[1];
					nmres[2] = nm[2];
					return getTriangleArea(triid);
				}
			
				return -1.0;
			}


		template<typename InputType>
		void  construct_from(BlackMesh<InputType>& pmesh)
		{
			clear();
		
			for (int i = 0; i < pmesh.GetNumVertices(); i++) {
				std::vector<Real> tmp;

				tmp.push_back((Real)Utils::to_double(pmesh.GetVertices()[i].pos[0]));
				tmp.push_back((Real)Utils::to_double(pmesh.GetVertices()[i].pos[1]));
				tmp.push_back((Real)Utils::to_double(pmesh.GetVertices()[i].pos[2]));
				insert_vtx(tmp);
			}

			auto itf = pmesh.GetTriangles();
			for (int i = 0; i < pmesh.GetNumTriangles(); i++) {
				auto it = itf[i].vertices;
				insert_face(it);
			}

			update_mesh_properties();
		}

		void update_mesh_to_keep_list(const std::vector<bool> &cc_should_keep)//;
		{
				auto tris = GetTrianglesCopy();
				clear_but_keep_vtx();
			
				for (int i = 0; i < tris.size(); i++) {
					if (cc_should_keep[tris[i].component_id]) {
						insert_face(tris[i].vertices);
					}
				}
			}
	};

	template <class Real>
	class Dual_graph : public BlackMesh<Real> {
	public: 
		class Edge{
		public:
			std::pair<int, int> vertices;//first should < second
			std::vector<int> incident_faces;
			bool is_manifold_edge = -1;
			int based_edge=-1;
		};

		inline void clear() {
			mNumVertices = mNumEdges = mNumTriangles = 0;
			mVertices.clear();
			mEdges.clear();
			mTriangles.clear();

			mAveEdge = 0;

			stacked_V.resize(0, 0);
			stacked_F.resize(0, 0);
			stacked_face_center.resize(0, 0);
			face_color_.resize(0, 0);
			stacked_E.resize(0, 0);
		}

		bool insert_edge(int i, int j, int dual_edge_id)//;//only should be used for dual graph where two nodes have multiple edges

		{
			//need to maintain the incident structure in Vertex
			assert(i != j);
			int query_id_min, query_id_max;
			if (i > j) { query_id_max = i; query_id_min = j; }
			else { query_id_max = j; query_id_min = i; }

			//query if this edge already exist
			int edge_exist = if_edge_already_exist(query_id_min, query_id_max, dual_edge_id);

			if (edge_exist != -1)
				return false;
			else {
				Dual_graph::Edge tmp;
				tmp.vertices.first = query_id_min;
				tmp.vertices.second = query_id_max;
				tmp.based_edge = dual_edge_id;

				mEdges.push_back(tmp);
				mNumEdges++;

				mVertices[i].incident_edges.push_back(mNumEdges - 1);
				mVertices[j].incident_edges.push_back(mNumEdges - 1);
			}


		}
		int if_edge_already_exist(int i, int j, int base_id)//;

		{
			assert(i != j);
			int query_id_min, query_id_max;
			if (i > j) { query_id_max = i; query_id_min = j; }
			else { query_id_max = j; query_id_min = i; }


			for (int k = 0; k < mVertices[i].incident_edges.size(); k++) {
				int edgeid = mVertices[i].incident_edges[k];
				Edge tmp = mEdges[edgeid];
				if (tmp.vertices.first == query_id_min&&tmp.vertices.second == query_id_max&&tmp.based_edge == base_id) {
					return edgeid;
				}

			}

			return -1;
		}

		std::vector<Edge> mEdges;
		inline const Edge* GetEdges() const { return &(mEdges[0]); }

		Eigen::MatrixXi stacked_E;
		void assemble_stacked_V_and_E()//;

		{
			stacked_V.resize(3, mNumVertices);
			for (int i = 0; i < mNumVertices; i++) {
				stacked_V.col(i) = Eigen::Vector3f(Utils::to_double(mVertices[i].pos[0]),
					Utils::to_double(mVertices[i].pos[1]), Utils::to_double(mVertices[i].pos[2]));
			}

			stacked_E.resize(2, mNumEdges);
			for (int i = 0; i < mNumEdges; i++) {
				stacked_E.col(i) = Eigen::Vector2i(mEdges[i].vertices.first, mEdges[i].vertices.second);
			}
		}
		void update_mesh_properties()//;
		{
			assert(mNumEdges = mEdges.size());
			assert(mNumVertices = mVertices.size());
			compute_avg_edge();
			assemble_stacked_V_and_E();
		}
	};


}


