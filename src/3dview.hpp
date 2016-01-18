#ifndef __3DVIEW_HPP__
#define __3DVIEW_HPP__

#define GL_GLEXT_PROTOTYPES 
#include <QtGui>

// define is needed for VBO function calls
#define GL_GLEXT_PROTOTYPES 
#include <QtOpenGL>
#include <QTreeWidget>
#include <QColor>

#include <vector>
#include <iostream>

#include "edge.hpp"

using namespace std;
using namespace geom;
using namespace vol;

class _3dView : public QGLWidget {
  Q_OBJECT
public:
  _3dView(QWidget *parent = 0);
  ~_3dView() {
    if(xz_img != NULL) { deleteTexture(xz_texture); delete xz_img; delete xz_fbo; }
    if(yz_img != NULL) { deleteTexture(yz_texture); delete yz_img; delete yz_fbo; }
    if(xy_img != NULL) { deleteTexture(xy_texture); delete xy_img; delete xy_fbo; }
  }
  void updateVol(edge::image_stack *newstack) {  
    image_stack = newstack;
    selected_part = selected_cell = NULL;

    setXYplane(xy_zpos); setYZplane(yz_xpos); setXZplane(xz_ypos);

    int vwidth = image_stack->width(), vheight = image_stack->height(), vdepth = image_stack->depth();
    center[0] = float(vwidth) / 2; center[1] = float(vheight) / 2; center[2] = float(vdepth) / 2;
  }
  void setVol(edge::image_stack *newstack) {  
    if(xz_img != NULL) { delete xz_img; xz_img = NULL; delete xz_fbo; xz_fbo = NULL; }
    if(yz_img != NULL) { delete yz_img; yz_img = NULL; delete yz_fbo; yz_fbo = NULL; }
    if(xy_img != NULL) { delete xy_img; xy_img = NULL; delete xy_fbo; xy_fbo = NULL; }
    
    image_stack = newstack;
    selected_cell = NULL;
    displayid = 0;

    int vwidth = image_stack->width(), vheight = image_stack->height(), vdepth = image_stack->depth();
    setXYplane(vdepth/2); setYZplane(vwidth/2); setXZplane(vheight/2);

    center[0] = float(vwidth) / 2; center[1] = float(vheight) / 2; center[2] = float(vdepth) / 2;

    // Position camera 3 * depth away from the image volume.
    camera.eye = center - 3 * vdepth * camera.towards;

    int maxdim = max(max(vwidth, vheight), vdepth);

    // NOTE: These have to be close to avoid artificats in depth buffer. 
    // TODO: Figure out how to set this correctly for given volume sizes.
    camera.neardist = 0.9 * (double)maxdim; camera.fardist = 5.5 * (double)maxdim;
  }
  void saveGLState() {
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glMatrixMode(GL_PROJECTION); glPushMatrix(); glMatrixMode(GL_MODELVIEW);  glPushMatrix();
  }
  void restoreGLState() {
    glMatrixMode(GL_PROJECTION); glPopMatrix();  glMatrixMode(GL_MODELVIEW);  glPopMatrix();
    glPopAttrib();
  }
  // Update position of XY plane
  void setXYplane(int newzpos) {
    xy_zpos = newzpos;
    vector<volume8 *> dispvols; getDispVols(dispvols);
    if(dispvols.size() == 0) return;
    volume8 *dispvol = dispvols[0];
    volume8 &v = *dispvol;
    if(xy_img == NULL) { // Nothing allocated yet.
      xy_img = new QImage(v.width, v.height, QImage::Format_ARGB32);
      xy_fbo = new QGLFramebufferObject(v.width, v.height);
    }
    else { // Handle resize.
      if(v.width != xy_img->width() || v.height != xy_img->height()) {
	delete xy_img; delete xy_fbo;
	xy_img = new QImage(v.width, v.height, QImage::Format_ARGB32);
	xy_fbo = new QGLFramebufferObject(v.width, v.height);
      }
      if(!image_stack->conf.use_GL_buffers) deleteTexture(xy_texture); // Use texture instead of FBO.
    }
    // Update image.
    if(dispvols.size() == 1) {
      for(int y = 0; y < v.height; y++) {
	QRgb *scanLine = (QRgb*)xy_img->scanLine(y);      
	for(int x = 0; x < v.width; x++) { unsigned char val = v(x, y, newzpos); scanLine[x] = qRgb(val,val,val);  }
      }
    }
    else if(dispvols.size() == 2) {
      volume8 &v2 = *dispvols[1];
      if(toggleColorsFlag) {
	for(int y = 0; y < v.height; y++) {
	  QRgb *scanLine = (QRgb*)xy_img->scanLine(y);      
	  for(int x = 0; x < v.width; x++) { 
	    unsigned char val2 = v(x, y, newzpos), val1 = v2(x, y, newzpos); 
	    scanLine[x] = qRgb(val2,val1,val2);  
	  }
	}
      }
      else {
	for(int y = 0; y < v.height; y++) {
	  QRgb *scanLine = (QRgb*)xy_img->scanLine(y);      
	  for(int x = 0; x < v.width; x++) { 
	    unsigned char val1 = v(x, y, newzpos), val2 = v2(x, y, newzpos); 
	    scanLine[x] = qRgb(val2,val1,val2);  
	  }
	}
      }
    }
    else if(dispvols.size() == 3) {
      volume8 &v2 = *dispvols[1], &v3 = *dispvols[2];
      for(int y = 0; y < v.height; y++) {
	QRgb *scanLine = (QRgb*)xy_img->scanLine(y);      
	for(int x = 0; x < v.width; x++) { 
	  unsigned char val1 = v(x, y, newzpos), val2 = v2(x, y, newzpos), val3 = v3(x, y, newzpos); 
	  scanLine[x] = qRgb(val2,val3,val1);  
	}
      }
    }    
    // Draw image on FBO.
    if(image_stack->conf.use_GL_buffers) {
      saveGLState();
      QPainter p(xy_fbo); p.drawImage(0,0, *xy_img); p.end();
      restoreGLState();
    }
    else xy_texture = bindTexture(*xy_img);
  }

  void setYZplane(int newxpos) {
    yz_xpos = newxpos;
    vector<volume8 *> dispvols; getDispVols(dispvols);
    if(dispvols.size() == 0) return;
    volume8 *dispvol = dispvols[0];
    volume8 &v = *dispvol;
    if(yz_img == NULL) {
      yz_img = new QImage(v.height, v.depth, QImage::Format_ARGB32);
      yz_fbo = new QGLFramebufferObject(v.height, v.depth);
    }
    else {
      if(yz_img->width() != v.height || yz_img->height() != v.depth) {
	delete yz_img; delete yz_fbo;
	yz_img = new QImage(v.height, v.depth, QImage::Format_ARGB32);
	yz_fbo = new QGLFramebufferObject(v.height, v.depth);
      }
      if(!image_stack->conf.use_GL_buffers) deleteTexture(yz_texture);
    }
    if(dispvols.size() == 1 ) {
      for(int z = 0; z < v.depth; z++) {
	QRgb *scanLine = (QRgb*)yz_img->scanLine(z);      
	for(int y = 0; y < v.height; y++) {
	  unsigned char val = v(newxpos, y, z);
	  scanLine[y] = qRgb(val,val,val);
	}
      }
    }
    else if(dispvols.size() == 2) {
      volume8 &v2 = *dispvols[1];
      if(toggleColorsFlag) {
	for(int z = 0; z < v.depth; z++) {
	  QRgb *scanLine = (QRgb*)yz_img->scanLine(z);      
	  for(int y = 0; y < v.height; y++) {
	    unsigned char val2 = v(newxpos, y, z), val1 = v2(newxpos, y, z);
	    scanLine[y] = qRgb(val2,val1,val2);
	  }
	}
      }
      else {
	for(int z = 0; z < v.depth; z++) {
	  QRgb *scanLine = (QRgb*)yz_img->scanLine(z);      
	  for(int y = 0; y < v.height; y++) {
	    unsigned char val1 = v(newxpos, y, z), val2 = v2(newxpos, y, z);
	    scanLine[y] = qRgb(val2,val1,val2);
	  }
	}
      }
    }
    else if(dispvols.size() == 3) {
      volume8 &v2 = *dispvols[1], &v3 = *dispvols[2];
      for(int z = 0; z < v.depth; z++) {
	QRgb *scanLine = (QRgb*)yz_img->scanLine(z);      
	for(int y = 0; y < v.height; y++) {
	  unsigned char val1 = v(newxpos, y, z), val2 = v2(newxpos, y, z), val3 = v3(newxpos, y, z);
	  scanLine[y] = qRgb(val2,val3,val1);
	}
      }
    }
    if(image_stack->conf.use_GL_buffers) {
      saveGLState();
      QPainter p(yz_fbo); p.drawImage(0,0, *yz_img); p.end();
      restoreGLState();
    }
    else yz_texture = bindTexture(*yz_img);
  }
  void setXZplane(int newypos) {
    xz_ypos = newypos;
    vector<volume8 *> dispvols; getDispVols(dispvols);
    if(dispvols.size() == 0) return;
    volume8 *dispvol = dispvols[0];
    volume8 &v = *dispvol;
    if(xz_img == NULL) {
      xz_img = new QImage(v.width, v.depth, QImage::Format_ARGB32);
      xz_fbo = new QGLFramebufferObject(v.width, v.depth);
    }
    else {
      if(xz_img->width() != v.width || xz_img->height() != v.depth) {
	delete xz_img; delete xz_fbo;
 	xz_img = new QImage(v.width, v.depth, QImage::Format_ARGB32);
	xz_fbo = new QGLFramebufferObject(v.width, v.depth);
      }
      if(!image_stack->conf.use_GL_buffers) deleteTexture(xz_texture);
    }
    if(dispvols.size() == 1) {
      for(int z = 0; z < v.depth; z++) {
	QRgb *scanLine = (QRgb*)xz_img->scanLine(z);      
	for(int x = 0; x < v.width; x++) {
	  unsigned char val = v(x, newypos, z);
	  scanLine[x] = qRgb(val,val,val);
	}
      }
    }
    else if(dispvols.size() == 2) {
      volume8 &v2 = *dispvols[1];
      if(toggleColorsFlag) {
	for(int z = 0; z < v.depth; z++) {
	  QRgb *scanLine = (QRgb*)xz_img->scanLine(z);      
	  for(int x = 0; x < v.width; x++) {
	    unsigned char val2 = v(x, newypos, z), val1 = v2(x, newypos, z);
	    scanLine[x] = qRgb(val2,val1,val2);
	  }
	}
      }
      else {
	for(int z = 0; z < v.depth; z++) {
	  QRgb *scanLine = (QRgb*)xz_img->scanLine(z);      
	  for(int x = 0; x < v.width; x++) {
	    unsigned char val1 = v(x, newypos, z), val2 = v2(x, newypos, z);
	    scanLine[x] = qRgb(val2,val1,val2);
	  }
	}
      }
    }
    else if(dispvols.size() == 3) {
      volume8 &v2 = *dispvols[1], &v3 = *dispvols[2];
      for(int z = 0; z < v.depth; z++) {
	QRgb *scanLine = (QRgb*)xz_img->scanLine(z);      
	for(int x = 0; x < v.width; x++) {
	  unsigned char val1 = v(x, newypos, z), val2 = v2(x, newypos, z), val3 = v3(x, newypos, z);
	  scanLine[x] = qRgb(val2,val3,val1);
	}
      }
    }
    if(image_stack->conf.use_GL_buffers) {
      saveGLState();
      QPainter p(xz_fbo); p.drawImage(0,0, *xz_img); p.end();
      restoreGLState();
    }
    else xz_texture = bindTexture(*xz_img);
  }

  // Returns the display volume based on current display ID.
  volume8 *getDispvol() { if(image_stack != NULL) return image_stack->dispvols[displayid]; else return NULL; }
  
  // If display ID is less than the number of channels.
  void getDispVols(vector<volume8 *> &vols) {
    if(image_stack == NULL) return;
    edgeconfig &conf = image_stack->conf;
    if(conf.keep_edge_channel_only || displayid >= image_stack->ichannels || image_stack->ichannels > 3  ) 
      vols.push_back(image_stack->dispvols[displayid]);
    else if(displayid < image_stack->ichannels) {
      for(int c = 0; c < image_stack->ichannels; c++) vols.push_back(image_stack->dispvols[c]);
    }
  }

  void toggleDOG() { // Cycle through intermediate processing volumes.
    if(image_stack == NULL) return;
    int tries = 0, nvols = image_stack->dispvols.size();
    while(true) {
      displayid++; if(displayid == nvols) displayid = 0;
      if(getDispvol() != NULL) break;
      if(tries == nvols) return;
      tries++;
    } 
    setYZplane(yz_xpos); setXYplane(xy_zpos); setXZplane(xz_ypos);
  }
  void toggleDOGrev() { // Cycle through intermediate processing volumes.
    if(image_stack == NULL) return;
    int tries = 0, nvols = image_stack->dispvols.size();
    while(true) {
      displayid--; if(displayid < 0) displayid = nvols-1;
      if(getDispvol() != NULL) break;
      if(tries == nvols) return;
      tries++;
    } 
    cout << "displayid = " << displayid << endl;     
    setYZplane(yz_xpos); setXYplane(xy_zpos); setXZplane(xz_ypos);
  }
  void toggleBgColor() { bgColorChange = true; bgColorWhite = !bgColorWhite; }

  // Sets viewer into rotate, scale, and translate mode. Works better on a track pad.
  void toggleRotate() { rotateFlag = !rotateFlag; setMouseTracking(rotateFlag); lastPosFlag = false; }
  void toggleScale() { scaleFlag = !scaleFlag; setMouseTracking(scaleFlag); lastPosFlag = false; }
  void toggleTranslate() { translateFlag = !translateFlag; setMouseTracking(translateFlag);  lastPosFlag = false; }
  void toggleUndisplay() { undisplayFlag = !undisplayFlag; }
  void toggleSelected() { selectFlag = !selectFlag; }
  void toggleTrajectoryMode() { 
    trajectoryMode = !trajectoryMode; 
    /*if(trajectoryMode && selected_cell != NULL) {
      quickselect_traj.clear();
      trajectoryid = selected_cell->trajectoryid;
      selected_cell->selected = false;
      selected_cell->displayed = true; 
      selected_cell = NULL; 
    }
    else if(trajectoryMode && quickselect_traj.size() > 0) trajectoryid = -1;*/
  }
  void toggleNeighborMode() {  neighborMode = !neighborMode; }
  void setQuickSelectTraj(vector<int> &selected) { quickselect_traj = selected; }
  void getQuickSelectTraj(vector<int> &out) { out = quickselect_traj; }
  bool isTrajectoryMode() { return trajectoryMode; }
  //int getTrajectoryId() { return trajectoryid; }
  void toggleEdit() { editFlag = !editFlag;  }
  bool getEditFlag() { return editFlag; }
  
  void toggleDisplay() { displayFlag = !displayFlag; }
  void toggleSlices() { slicesFlag = !slicesFlag; }
  void togglePolygonMode() { polygonMode = !polygonMode; }
  void toggleColors() { 
    toggleColorsFlag = !toggleColorsFlag;
    if(image_stack == NULL) return;
    setYZplane(yz_xpos); setXYplane(xy_zpos); setXZplane(xz_ypos);
  }
  void toggleDisplayNucs() { displayNucsFlag = !displayNucsFlag; }

  void clearSelected() { if(selected_cell != NULL) selected_cell->displayed = true; selected_cell = NULL; }

  edge::image_cell *getSelectedCell() { return selected_cell; }
  void setMeasurementTree(QTreeWidget *zTree) { this->zTree = zTree; }

  void toggleSamples() { if(selected_cell != NULL) displaySamples = !displaySamples; }
  void toggleLabels() { displayLabels = !displayLabels; }
  void cycleLabels() { if(selected_cell != NULL) selected_cell->cycle_labels();  }

  void updateMeasurementTree() { if(selected_cell != NULL) updateMeasurementTree(selected_cell); }

  void setProject(edge::hyperstack *p) { project = p; }
  void getPos(double &yz_x, double &xz_y, double &xy_z) { yz_x = yz_xpos; xz_y = xz_ypos; xy_z = xy_zpos;  }
  void toggleBfsMode() { bfsMode = !bfsMode; }
  void add_selected_to_quickselect() {
    if(selected_cell == NULL) { cout << "no selected cell" << endl; return; }
    int sel_traj_id = selected_cell->trajectoryid;
    cout << "add sel = " << sel_traj_id << endl;
    bool found = false;
    for(size_t q = 0; q < quickselect_traj.size(); q++) {
      if(quickselect_traj[q] == sel_traj_id) { found = true; break; }
    }
    if(!found) quickselect_traj.push_back(sel_traj_id);
  }
  void remove_selected_from_quickselect() {
    if(selected_cell == NULL) return;
    int sel_traj_id = selected_cell->trajectoryid;
    for(size_t q = 0; q < quickselect_traj.size(); q++) {
      if(quickselect_traj[q] == sel_traj_id) {
	quickselect_traj[q] = quickselect_traj.back();
	quickselect_traj.pop_back();
	return;
      }
    }
  }
protected:
  void initializeGL(); // QGLWidget OpenGL interface
  void paintGL();
  void resizeGL(int width, int height);

  // OpenGL init and draw routines.
  void initCameraGL();
  void drawStackGL();
  void drawSelectedGL();
  void drawCellsGL();
  void updateMeasurementTree(edge::image_cell *cell);
  void drawCellGL(edge::image_cell *cell, GLfloat front[4], GLfloat back[4], bool skip_gl_buffers = false);
  void drawCellLabeled(edge::image_cell *cell);
  void drawCellLabeledNeighbors(edge::image_cell *cell, vector<edge::image_cell *> &cells);

  void addItem2Tree(const QString key, QString value, QString units);
  void updateItem(int itemNum, QString value);

  void mousePressEvent(QMouseEvent *event); void mouseMoveEvent(QMouseEvent *event);
private:
  edge::hyperstack *project;
  int cur_stack;

  QImage *xy_img, *yz_img, *xz_img;  // Corresponding QImage objects for volume slices.
  QGLFramebufferObject *xy_fbo, *yz_fbo, *xz_fbo; // FBOs for displaying slices.

  double xy_zpos, yz_xpos, xz_ypos;
  int displayid; //trajectoryid;
  vector<int> quickselect_traj;

  // Flags that encode GUI state information.
  bool bfsMode;
  bool editFlag, displayFlag;
  bool slicesFlag;
  bool selectFlag, undisplayFlag, scaleFlag, translateFlag, rotateFlag;
  bool bgColorChange, bgColorWhite;
  bool displayLabels, displaySamples;
  bool trajectoryMode, neighborMode;
  bool toggleColorsFlag;
  bool displayNucsFlag;
  bool polygonMode;

  // 3dview geometry 
  camera3 camera; vec3 center;

  // Used to track mouse position for keyboard shortcuts.
  bool lastPosFlag; // Flag is true if a last position has been captured. Set to false 
                    // prior to mouse driven interaction.
  QPoint lastPos;   // Mouse state information.

  GLuint xy_texture, yz_texture, xz_texture;

  edge::image_stack *image_stack;
  edge::image_cell *selected_cell, *selected_part; 
  QTreeWidget *zTree;

  vector<QColor> colormap, colormap_bfs;
};

#endif

