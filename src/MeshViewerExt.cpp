#include "stdafx.h"
#include "MeshViewerExt.h"

Optimizer opt;

MeshViewerExt::MeshViewerExt(const char* _title, int _width, int _height)
	: MeshViewer(_title, _width, _height), show_selection_(true), is_kdTree_init_(false) {
	obj_sphere_ = gluNewQuadric();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	gluQuadricDrawStyle(obj_sphere_, GLU_FILL);
	gluQuadricNormals(obj_sphere_, GLU_SMOOTH);

	m_flood_type_ = BY_MANIFOLD;
	opt_ = &opt;

	min_timer_ = 1;
	max_timer_ = 15;

	sWeight_ = 0.1;
	uWeight_ = 1.0;
	bWeight_ = 1.0;
}

MeshViewerExt::~MeshViewerExt() {
	gluDeleteQuadric(obj_sphere_);
}


