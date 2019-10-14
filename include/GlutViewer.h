#pragma once

//#include <AntTweakBar.h>
#include <gl/freeglut.h>

class GlutViewer {
 public:
  enum DrawMode { HIDDEN_LINE = 1, WIRE_FRAME, SOLID_FLAT, SOLID_SMOOTH };
  enum ThemeType { DARK_THEME = 1, LIGHT_THEME };
  enum MaterialType { BRASS = 1, COPPER, GOLD, JADE, LIGHT_BLUE, SILVER, PEWTER, OBSIDIAN, PURPLE, TURQUOISE };

 public:
  GlutViewer(const char* _title, int _width, int _height);
  virtual ~GlutViewer();
  void launch(void);

 protected:
  virtual void setup_glut(void);
  virtual void setup_view(void);
  virtual void setup_scene(const Float* _cog, const Float _radius);

  void apply_projection_matrix(bool draw_3d = true);
  void apply_modelview_matrix(bool draw_3d = true);



 protected:
  // screen width and height and title
  int  width_, height_;
  string title_;

  // scene position and dimension
  double center_[3];
  double radius_;
  double trans_[3];
  double rotation_matrix_[16];
  double scale_;

  // projection parameters
  double near_, far_, fovy_;

  // OpenGL matrices
  GLdouble projection_matrix_[16];
  GLdouble modelview_matrix_[16];
  GLint viewport_[4];

  // trackball helpers
  int last_P2d_[2];
  double last_P3d_[3];
  bool button_down_[5];

  // draw mode
  DrawMode draw_mode_;

  // Theme type
  ThemeType theme_type_;
  float bak_color_[3];
  float line_color_[3];

  // material type
  MaterialType material_type_;

  // light direction
  float light_direction_[4];

  // light model
  bool two_side_rendering_;
  bool local_viewer_;


  // for full screen
  bool fullscreen_;
  int  bak_left_, bak_top_, bak_width_, bak_height_;

  // for timer
  unsigned int msecs_;

 private:
  static GlutViewer* current_viewer_;
};
