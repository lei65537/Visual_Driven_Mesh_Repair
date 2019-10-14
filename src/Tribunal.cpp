#include "stdafx.h"
#include "Tribunal.h"



void Tribunal::assemble_number_of_edge()
{
	if (_gh->non_manifold_curves.size() == 0) {
		std::cout << "Curve not assembled yet\n";
#if defined(_RECORD_ERROR_LOGS) 
		(*_log_off) << _file_under_processing << " ERROR Curve not assembled yet\n";
		(*_log_off).flush();
#endif
	}

	cc_incident_edge_num.resize(_gh->node_list.size(), 0);
	for (int i = 0; i < _gh->non_manifold_curves.size(); i++) {
		int someedgeid = _gh->non_manifold_curves[i][0];

		//count number
		for (int j = 0; j < _gh->edge_list[someedgeid].incident_components.size(); j++) {
			int ccid = _gh->edge_list[someedgeid].incident_components[j];
			cc_incident_edge_num[ccid]++;
		}
	}
}

Tribunal::Tribunal()
{
}


Tribunal::~Tribunal()
{
}




void Tribunal::build_cgal_mesh() {
	_triangles.clear();
	_triangles.reserve(3 * _msh->GetNumTriangles());

	for (int i = 0; i < _msh->GetNumTriangles(); i++) {
		_triangles.push_back(My_point(_msh->GetVertices()[_msh->GetTriangles()[i].vertices[0]].pos));
		_triangles.push_back(My_point(_msh->GetVertices()[_msh->GetTriangles()[i].vertices[1]].pos));
		_triangles.push_back(My_point(_msh->GetVertices()[_msh->GetTriangles()[i].vertices[2]].pos));
	}
}

void Tribunal::build_AABB_Tree() {
	_tree.insert(Triangle_iterator(_triangles.begin()),
		Triangle_iterator(_triangles.end()));
	_tree.build();
}

bool Tribunal::if_intersects(array<double, 3>& pos, std::vector<double> &nm, double CUTTING_EPS)
{
	//intersect - return true
	K::Ray_3 ray_query(K::Point_3(pos[0] + CUTTING_EPS*nm[0], pos[1] + CUTTING_EPS*nm[1], pos[2] + CUTTING_EPS*nm[2]),
		K::Point_3(pos[0] + nm[0], pos[1] + nm[1], pos[2] + nm[2]));
	return _tree.number_of_intersected_primitives(ray_query) != 0;
}


void Tribunal::update_unary_score(double CUTTING_EPS, double IGNOR_TRI_UNDER,
	std::vector<double>* pt_posSco, std::vector<double>* pt_negSco)
{

	std::vector<int> posPix, negPix;
	VisualProcesser vp;
	vp.setMesh(_msh);
	vp.get_unary_score_all_campos(posPix, negPix);

	std::vector<int> posPixCC, negPixCC;
	posPixCC.resize(_sampler->cc_triid.size(), 0);
	negPixCC.resize(_sampler->cc_triid.size(), 0);

	for (int i = 0; i < _msh->GetNumTriangles(); i++) {
		int ccid = _msh->GetTriangles()[i].component_id;
		posPixCC[ccid] += posPix[i];
		negPixCC[ccid] += negPix[i];
	}

	for (int i = 0; i < _sampler->cc_triid.size(); i++) {
		int totalPix = posPixCC[i] + negPixCC[i];
		if (totalPix > _TREAT_PIX_SMALLER_THAN_AS_ZERO) {
			double score_pos = (double)posPixCC[i] / (double)totalPix;
			double score_neg = (double)negPixCC[i] / (double)totalPix;

			if (_unary_score != NULL) {
				if (cc_incident_edge_num[i] > 0) {
					(*_unary_score)[i].pos_score = score_pos / (double)cc_incident_edge_num[i];
					(*_unary_score)[i].neg_score = score_neg / (double)cc_incident_edge_num[i];
				}
				else {
					(*_unary_score)[i].pos_score = score_pos ;
					(*_unary_score)[i].neg_score = score_neg ;
				}
			}

			if (pt_posSco != NULL&&pt_negSco != NULL) {
				(*pt_posSco)[i] = score_pos;
				(*pt_negSco)[i] = score_neg;
			}
		}
		else {
			if (_unary_score != NULL) {
				if (cc_incident_edge_num[i] > 0) {
					(*_unary_score)[i].pos_score = 0.5 / (double)cc_incident_edge_num[i];
					(*_unary_score)[i].neg_score = 0.5 / (double)cc_incident_edge_num[i];
				}
				else {
					(*_unary_score)[i].pos_score = 0.5 ;
					(*_unary_score)[i].neg_score = 0.5 ;
				}
			}
			if (pt_posSco != NULL&&pt_negSco != NULL) {
				(*pt_posSco)[i] = 0.5;
				(*pt_negSco)[i] = 0.5;
			}
		}
	}


}

void Tribunal::upate_binary_score(double IGNOR_TRI_UNDER) {
	//calculate length
	std::vector<double> non_manifold_edge_len;
	non_manifold_edge_len.resize(_gh->edge_list.size(), 0.0);
	std::vector<double> curve_len;
	curve_len.resize(_gh->non_manifold_curves.size(), 0.0);

	for (int i = 0; i < _gh->non_manifold_curves.size(); i++) {
		for (int j = 0; j < _gh->non_manifold_curves[i].size(); j++) {
			int edgeid = _gh->non_manifold_curves[i][j];
			auto ehgh = _gh->edge_list[edgeid];
			auto ehmsh = _msh->GetEdges()[ehgh.original_edge_id];
			auto vtxs = _msh->GetVertices()[ehmsh.vertices.first];
			auto vtxe = _msh->GetVertices()[ehmsh.vertices.second];

			double len = edge_length(&(vtxs.pos[0]), &(vtxe.pos[0]));
			non_manifold_edge_len[edgeid] = len;
			curve_len[i] += len;
		}
	}

	for (auto it = _binary_score->begin(); it != _binary_score->end(); it++) {
		int curvid = it->first.first;
		//int dualedgeid = it->first.second;
		int indi = it->second.inidi;
		int indj = it->second.inidj;
		int dirfi = it->second.diri;
		int dirfj = it->second.dirj;

#ifdef _PUNISH_CONNECTION_ON_CURV_SMALLER_THAN
		bool punish_flag = false;
		{
			auto eh = _gh->non_manifold_curves[curvid][0];
			int ccidi = _gh->edge_list[eh].incident_components[indi];
			int ccidj = _gh->edge_list[eh].incident_components[indj];
			double cv_to_area = curve_len[curvid] / 
				sqrt(0.5*(_sampler->cc_total_area[ccidi] + _sampler->cc_total_area[ccidj]));
			if (cv_to_area < _PUNISH_CONNECTION_ON_CURV_SMALLER_THAN) {
				punish_flag = true;
			}
		}

		if (!punish_flag) {
#endif
			double weighted_sum = 0.0;
			double acc_len = 0.0;
			for (int j = 0; j < _gh->non_manifold_curves[curvid].size(); j++) {
				//all edges in this curve
				int edgeid = _gh->non_manifold_curves[curvid][j];
				auto ehgh = _gh->edge_list[edgeid];
				auto ehmsh = _msh->GetEdges()[ehgh.original_edge_id];

				int idfi = ehmsh.incident_faces[indi];
				int idfj = ehmsh.incident_faces[indj];


				if (_sampler->face_area[idfi] < IGNOR_TRI_UNDER ||
					_sampler->face_area[idfj] < IGNOR_TRI_UNDER) {


					double hei1[3];
					double hei2[3];
					double flag1 = _msh->compute_fake_normal(idfi, ehgh.original_edge_id, hei1, 0, _SEARCH_NEIB_ON_DEGENERATE);
					double flag2 = _msh->compute_fake_normal(idfj, ehgh.original_edge_id, hei2, 0, _SEARCH_NEIB_ON_DEGENERATE);

					if (flag1 < 0 || flag2 < 0) {
						std::cout << "DEG_TRI_FOUND\n";
#if defined(_RECORD_ERROR_LOGS) 
						(*_log_off) << _file_under_processing << " WARNING degenerate found in binary term\n";
						(*_log_off).flush();
#endif
						continue;
					}

					//hei1 and hei2 already normalized

					double scli = dirfi == 0 ? 1.0 : -1.0;
					double sclj = dirfj == 0 ? 1.0 : -1.0;

					auto tmp1 = std::vector<double>{ scli*hei1[0],scli*hei1[1] ,scli*hei1[2] };
					auto tmp2 = std::vector<double>{ sclj*hei2[0],sclj*hei2[1] ,sclj*hei2[2] };

					double cos_sco = (cos_angle(&(tmp1[0]), &(tmp2[0])));

					if (cos_sco > 1.0 - 1e-10) {
						cos_sco = 1.0 - 1e-10;
					}
					else if (cos_sco < -1.0 + 1e-10) {
						cos_sco = -1.0 + 1e-10;
					}
					cos_sco = 1.0 - acos(cos_sco) / M_PI;

					weighted_sum += non_manifold_edge_len[edgeid] * cos_sco;
					acc_len += non_manifold_edge_len[edgeid];

				}
				else {

					auto nm1 = _msh->GetTriangles()[idfi].normal;
					auto nm2 = _msh->GetTriangles()[idfj].normal;

					double scli = dirfi == 0 ? 1.0 : -1.0;
					double sclj = dirfj == 0 ? 1.0 : -1.0;

					auto tmp1 = std::vector<double>{ scli*nm1[0],scli*nm1[1] ,scli*nm1[2] };
					auto tmp2 = std::vector<double>{ sclj*nm2[0],sclj*nm2[1] ,sclj*nm2[2] };



					double cos_sco = (cos_angle(&(tmp1[0]), &(tmp2[0])));
					if (cos_sco > 1.0 - 1e-10) {
						cos_sco = 1.0 - 1e-10;
					}
					else if (cos_sco < -1.0 + 1e-10) {
						cos_sco = -1.0 + 1e-10;
					}
					cos_sco = 1.0 - acos(cos_sco) / M_PI;



					weighted_sum += non_manifold_edge_len[edgeid] * cos_sco;
					acc_len += non_manifold_edge_len[edgeid];
				}
			}

			if (acc_len > 0)
				it->second.score = weighted_sum / (acc_len);
			else
				it->second.score = 0;

#ifdef _PUNISH_CONNECTION_ON_CURV_SMALLER_THAN
		}
		else {
			it->second.score = _PUNISH_SCORE;
		}
#endif
	}


	
	for (auto it = _unary_score_from_binary->begin(); it != _unary_score_from_binary->end(); it++) {
		int curvid = it->first.first;
		//int dualedgeid = it->first.second;
		int indi = it->second.inidi;
		int indj = it->second.inidj;
		int dirfi = it->second.diri;
		//int dirfj = it->second.dirj;


		double weighted_sum = 0.0;
		double acc_len = 0.0;


			for (int j = 0; j < _gh->non_manifold_curves[curvid].size(); j++) {
				//all edges in this curve
				int edgeid = _gh->non_manifold_curves[curvid][j];
				auto ehgh = _gh->edge_list[edgeid];
				auto ehmsh = _msh->GetEdges()[ehgh.original_edge_id];

				int idfi = ehmsh.incident_faces[indi];
				int idfj = ehmsh.incident_faces[indj];

				if (_sampler->face_area[idfi] < IGNOR_TRI_UNDER ||
					_sampler->face_area[idfj] < IGNOR_TRI_UNDER) {

					double hei1[3];
					double hei2[3];
					double flag1 = _msh->compute_fake_normal(idfi, ehgh.original_edge_id, hei1, 0, _SEARCH_NEIB_ON_DEGENERATE);
					double flag2 = _msh->compute_fake_normal(idfj, ehgh.original_edge_id, hei2, 0, _SEARCH_NEIB_ON_DEGENERATE);

					if (flag1 < 0 || flag2 < 0) {
						std::cout << "DEG_TRI_FOUND\n";
#if defined(_RECORD_ERROR_LOGS) 
						(*_log_off) << _file_under_processing << " WARNING degenerate found in unary_binary \n";
						(*_log_off).flush();
#endif
						continue;
					}

					//hei1 and hei2 already normalized
					double scli = dirfi == 0 ? 1.0 : -1.0;
					//double sclj = dirfj == 0 ? 1.0 : -1.0;

					auto tmp1 = std::vector<double>{ scli*hei1[0],scli*hei1[1] ,scli*hei1[2] };
					auto tmp2 = std::vector<double>{ scli*hei2[0],scli*hei2[1] ,scli*hei2[2] };
					double cos_sco = (cos_angle(&(tmp1[0]), &(tmp2[0])));

					if (cos_sco > 1.0 - 1e-10) {
						cos_sco = 1.0 - 1e-10;
					}
					else if (cos_sco < -1.0 + 1e-10) {
						cos_sco = -1.0 + 1e-10;
					}
					cos_sco = 1.0 - acos(cos_sco) / M_PI;

					weighted_sum += non_manifold_edge_len[edgeid] * cos_sco;
					acc_len += non_manifold_edge_len[edgeid];

				}
				else {
					auto nm1 = _msh->GetTriangles()[idfi].normal;
					auto nm2 = _msh->GetTriangles()[idfj].normal;

					double scli = dirfi == 0 ? 1.0 : -1.0;
					//double sclj = dirfj == 0 ? 1.0 : -1.0;

					auto tmp1 = std::vector<double>{ scli*nm1[0],scli*nm1[1] ,scli*nm1[2] };
					auto tmp2 = std::vector<double>{ scli*nm2[0],scli*nm2[1] ,scli*nm2[2] };


					double cos_sco = (cos_angle(&(tmp1[0]), &(tmp2[0])));
					if (cos_sco > 1.0 - 1e-10) {
						cos_sco = 1.0 - 1e-10;
					}
					else if (cos_sco < -1.0 + 1e-10) {
						cos_sco = -1.0 + 1e-10;
					}
					cos_sco = 1.0 - acos(cos_sco) / M_PI;



					weighted_sum += non_manifold_edge_len[edgeid] * cos_sco;
					acc_len += non_manifold_edge_len[edgeid];
				}
			}

			if (acc_len > 0)
				it->second.score = weighted_sum / (acc_len);
			else
				it->second.score = 0;



	}
}

void Tribunal::update_weight(std::vector<double> &curve_len, double uw, double bw)
{

	//unary weight
	for (int i = 0; i < _unary_score->size(); i++) {
		(*_unary_score)[i].weight = uw*_sampler->cc_total_area[i] / _sampler->total_area;
	}

	total_dual_edge_weight = 0;

	

		double total_weight = 0.0;
		for (auto it = _binary_score->begin(); it != _binary_score->end(); it++) {
			int curveid = it->first.first;
			int ccidi = _gh->edge_list[_gh->non_manifold_curves[curveid][0]]
				.incident_components[it->second.inidi];
			int ccidj = _gh->edge_list[_gh->non_manifold_curves[curveid][0]]
				.incident_components[it->second.inidj];
			it->second.weight = curve_len[curveid];
			total_weight += it->second.weight;
		}

		//unary-binary weight
		for (auto it = _unary_score_from_binary->begin(); it != _unary_score_from_binary->end(); it++) {
			int curveid = it->first.first;
			int ccidi = _gh->edge_list[_gh->non_manifold_curves[curveid][0]]
				.incident_components[it->second.inidi];

			it->second.weight = curve_len[curveid];
			total_weight += it->second.weight;
		}


		//divided by total weight
		for (auto it = _binary_score->begin(); it != _binary_score->end(); it++) {
			it->second.weight = bw*it->second.weight / total_weight;
		}

		//unary-binary weight
		for (auto it = _unary_score_from_binary->begin(); it != _unary_score_from_binary->end(); it++) {
			it->second.weight = bw*it->second.weight / total_weight;
		}
		total_dual_edge_weight = total_weight;


}


void Tribunal::GroundColorMix(const double x_in, double & r, double & g, double & b)
{
	const double rone = 0.8;
	const double gone = 1.0;
	const double bone = 1.0;
	double x = x_in;
	x = (x_in < 0 ? 0 : (x > 1 ? 1 : x));

	if (x < 1. / 8.)
	{
		r = 0;
		g = 0;
		b = bone*(0.5 + (x) / (1. / 8.)*0.5);
	}
	else if (x < 3. / 8.)
	{
		r = 0;
		g = gone*(x - 1. / 8.) / (3. / 8. - 1. / 8.);
		b = bone;
	}
	else if (x < 5. / 8.)
	{
		r = rone*(x - 3. / 8.) / (5. / 8. - 3. / 8.);
		g = gone;
		b = (bone - (x - 3. / 8.) / (5. / 8. - 3. / 8.));
	}
	else if (x < 7. / 8.)
	{
		r = rone;
		g = (gone - (x - 5. / 8.) / (7. / 8. - 5. / 8.));
		b = 0;
	}
	else
	{
		r = (rone - (x - 7. / 8.) / (1. - 7. / 8.)*0.5);
		g = 0;
		b = 0;
	}
}

void Tribunal::visualize_unary_score(const char * filename)
{
	//here the score is  score*incident_edge_number

	char buffer1[128];
	sprintf(buffer1, "%s_pos.off", filename);
	ofstream offpos(buffer1);

	char buffer2[128];
	sprintf(buffer2, "%s_neg.off", filename);
	ofstream offneg(buffer2);

	offpos << "OFF\n";
	offpos << _msh->GetNumVertices() << " " << _msh->GetNumTriangles() << " 0\n";

	offneg << "OFF\n";
	offneg << _msh->GetNumVertices() << " " << _msh->GetNumTriangles() << " 0\n";
	for (int i = 0; i < _msh->GetNumVertices(); i++) {
		auto vh = _msh->GetVertices()[i];
		offpos << vh.pos[0] << " " << vh.pos[1] << " " << vh.pos[2] << " \n";
		offneg << vh.pos[0] << " " << vh.pos[1] << " " << vh.pos[2] << " \n";
	}

	for (int i = 0; i < _msh->GetNumTriangles(); i++) {
		auto fh = _msh->GetTriangles()[i];
		offpos << "3 " << fh.vertices[0] << " " << fh.vertices[1] << " " << fh.vertices[2] << " ";
		offneg << "3 " << fh.vertices[0] << " " << fh.vertices[1] << " " << fh.vertices[2] << " ";

		auto scpos = (*_unary_score)[fh.component_id].pos_score*cc_incident_edge_num[fh.component_id];
		double colorpos[3];
		GroundColorMix(scpos, colorpos[0], colorpos[1], colorpos[2]);
		offpos << (int)(colorpos[0] * 255) << " " << (int)(colorpos[1] * 255) << " " << (int)(colorpos[2] * 255) << " \n";

		auto scneg = (*_unary_score)[fh.component_id].neg_score*cc_incident_edge_num[fh.component_id];
		GroundColorMix(scneg, colorpos[0], colorpos[1], colorpos[2]);
		offneg << (int)(colorpos[0] * 255) << " " << (int)(colorpos[1] * 255) << " " << (int)(colorpos[2] * 255) << " \n";
	}

	offpos.close();
	offneg.close();
}
