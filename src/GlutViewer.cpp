#include "stdafx.h"
#include "GlutViewer.h"

// -----------
// static data component
GlutViewer* GlutViewer::current_viewer_ = NULL;

// -----------
// constructor and destructor
GlutViewer::GlutViewer(const char* _title, int _width, int _height)
  : title_(_title), width_(_width), height_(_height) {
  // screen center and radius
  scale_ = 1.0;
  radius_ = 1;
  for (int i = 0; i < 3; i++) {
    trans_[i] = center_[i] = 0.0;
  }

  // projection parameter
  near_ = 0.001 * radius_;
  far_ = 10.0 * radius_;
  fovy_ = 45.0;

  for (int i = 0; i < 16; i++) {
    double v = (i % 5) ? 0.0 : 1.0;
    rotation_matrix_[i] = v;
    projection_matrix_[i] = v;
    modelview_matrix_[i] = v;
  }

  // light direction
  light_direction_[0] = 0;
  light_direction_[1] = 0;
  light_direction_[2] = 1;
  light_direction_[3] = 0;

  // light model
  two_side_rendering_ = false;
  local_viewer_ = false;

  // menu related variables
  draw_mode_ = SOLID_FLAT;

  // material type
  //material_type_ = LIGHT_BLUE;
  material_type_ = PEWTER;

  // theme mode
  //theme_type_ = LIGHT_THEME;
  //bak_color_[0] = bak_color_[1] = bak_color_[2] = 1.0f;
  //line_color_[0] = line_color_[1] = line_color_[2] = 0.5f;
  theme_type_ = DARK_THEME;
  bak_color_[0] = bak_color_[1] = bak_color_[2] = 0.196f;
  line_color_[0] = line_color_[1] = line_color_[2] = 1.0f;

  // init mouse buttons
  for (int i = 0; i < 5; ++i)
    button_down_[i] = false;

  // not full screen
  fullscreen_ = false;

  // for timer
  msecs_ = 300;

  // pointer to current viewer
  current_viewer_ = this;
}

GlutViewer::~GlutViewer() {
  glutDestroyWindow(glutGetWindow());
}

// -----------
// public interface
void GlutViewer::launch(void) {
  setup_glut();
  setup_view();
  Float center[] = { 0.0, 0.0, 0.0 };
  setup_scene(center, 1);

}

// -----------
// setup
void GlutViewer::setup_glut(void) {
  // create window
  glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE | GLUT_ALPHA | GLUT_MULTISAMPLE);
  glutInitWindowSize(width_, height_);
  glutCreateWindow(title_.c_str());

}


void GlutViewer::setup_view(void) {
  // OpenGL state

  glClearColor(bak_color_[0], bak_color_[1], bak_color_[2], 0.0);
  //glDisable(GL_DITHER);
  glEnable(GL_DEPTH_TEST);

  // lighting
  glEnable(GL_LIGHT0);

  // anti-alias
  glEnable(GL_POINT_SMOOTH);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_POLYGON_SMOOTH);
  glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

  glEnable(GLUT_MULTISAMPLE);

  // blend
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}


void GlutViewer::setup_scene(const Float* _cog, const Float _radius) {
  for (int i = 0; i < 3; i++) {
    // since radius changed, the translation needs rescale here
    trans_[i] = trans_[i] * _radius / radius_;
    center_[i] = _cog[i];
  }

  radius_ = _radius;
  near_ = 0.001 * radius_;
  far_ = 10.0 * radius_;



  apply_projection_matrix();
  apply_modelview_matrix();
}

void GlutViewer::apply_projection_matrix(bool draw_3d) {
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  if (draw_3d) {
    gluPerspective(fovy_, (GLdouble)width_ / (GLdouble)height_, near_, far_);
    glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix_);
  } else {
    gluOrtho2D(0, width_, 0, height_);
  }
  glMatrixMode(GL_MODELVIEW);
}

void GlutViewer::apply_modelview_matrix(bool draw_3d) {
  // scene pos and size
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  if (draw_3d) {
    glTranslated(0, 0, -3 * radius_);					// view all
    glTranslated(trans_[0], trans_[1], trans_[2]);		// translate
    glMultMatrixd(rotation_matrix_);					// rotation
    //glScaled(scale_, scale_, scale_);					// scale
    glTranslated(-center_[0], -center_[1], -center_[2]);//	move the object to (0,0,0)
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix_);
  }
}
