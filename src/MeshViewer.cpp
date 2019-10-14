#include "stdafx.h"
#include "MeshViewer.h"
//#include <LibGP/write_png.h> 
#include <time.h>

#define EIGEN_MESHVIEWER
BlackMesh::BlackMesh<double> mesh;
BlackMesh::Dual_graph<double> dual_gh;
Preprocess<double> precesser;

MeshViewer::MeshViewer(const char* _title, int _width, int _height)
  : GlutViewer(_title, _width, _height), mesh_(nullptr), show_clip_box_(false), fix_clip_box_(false), update_bbox_(true),
    bbMin_(-0.7, -0.7, -0.7), bbMax_(0.7, 0.7, 0.7), use_face_color_(false), batch_screen_shot_(false), show_dual_graph_(false)
{
	show_mesh_ = true;
	show_dual_graph_ = false;

	//bar_rmsh_ = TwNewBar("Remeshing");
}

MeshViewer::~MeshViewer() {}
