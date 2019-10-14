#pragma once
#include "GlutViewer.h"
#include "BlackMesh.h"
#include "MeshIO.h"
#include "Preprocess.h"
//#include "TriMesh.h"

#if defined SINGLE_PRECISION
#define DEFINED_GL_FLOAT GL_FLOAT
#define  glNormal3Fv glNormal3fv
#define  glVertex3Fv glVertex3fv
#else
#define DEFINED_GL_FLOAT GL_DOUBLE
#define  glNormal3Fv glNormal3dv
#define glVertex3Fv glVertex3dv
#endif



enum FileType { MODEL_FILE = 1, IMAGE_FILE, COLOR_FILE };

class MeshViewer : public GlutViewer {
 public:
  // default constructor
  MeshViewer(const char* _title, int _width, int _height);
  ~MeshViewer();


 protected:
  // filenames
  string filename_; // full name
  string model_name_;
  string model_info_;
  string message_;
  vector<string> all_filenames_;

  // pointer of mesh
  BlackMesh::BlackMesh<double>* mesh_;
  BlackMesh::Dual_graph<double>* dual_gh_;
  Preprocess<double> *precesser_;

  // bounding points - for mesh drawing
  Vector3F bbMin_, bbMax_;

  // for face color
  bool use_face_color_;

  // for screen shot
  bool show_clip_box_;
  bool fix_clip_box_;
  bool update_bbox_;
  bool batch_screen_shot_;
  // screen shot bounding box
  Vector2i cbMin_, cbMax_;

  bool show_mesh_;
  bool show_dual_graph_;
};

