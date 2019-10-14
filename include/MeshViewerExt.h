#pragma once
#include "MeshViewer.h"
//#include <LibGP/KDTreeAdaptor.h>

#include "Optimizer.h"
#include "OffScreenRenderer.h"
#include "VisualProcesser.h"
#include "Intersector.h"
#include "Tribunal.h"

class  MeshViewerExt : public MeshViewer {
 public:
  // default constructor
  MeshViewerExt(const char* _title, int _width, int _height);
  ~MeshViewerExt();

 

 private:
  // for selecting point
 // LibGP::KDTree kdtree_pts_;
  //LibGP::KDTree kdtree_faces_;
  vector<int> selected_pts_;
  vector<int> selected_faces_;
  bool show_selection_;
  bool is_kdTree_init_;

  //
  Optimizer *opt_;

  // sphere
  GLUquadricObj* obj_sphere_;

  enum FloodType { BY_CONNECTION=1, BY_MANIFOLD };
  FloodType m_flood_type_;

  //
  int min_timer_;
  int max_timer_;

  //
  double sWeight_;
  double uWeight_;
  double bWeight_;
};