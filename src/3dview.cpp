#define GL_GLEXT_PROTOTYPES
#include<QtGui>

#include <algorithm>

#include "3dview.hpp"

using namespace std;

GLfloat grayGL[4] = {0.5, 0.5, 0.5, 1}; GLfloat greenGL[4] = {0, 1, 0, 1};
GLfloat blueGL[4] = {0, 0, 1, 1};       GLfloat redGL[4] = {1, 0, 0, 1};
GLfloat yellowGL[4] = {1, 1, 0, 1};     GLfloat magentaGL[4] = {1, 0, 1, 1};
GLfloat cyanGL[4] = {0, 1, 1, 1};

_3dView::_3dView(QWidget *parent)  : QGLWidget(parent) {
  project = NULL;
  image_stack = NULL;   // No image stack initially.
  selected_cell = selected_part = NULL; // No selected cell.

  xz_img = yz_img = xy_img = NULL;
  xz_fbo = yz_fbo = xy_fbo = NULL;

  editFlag = undisplayFlag = translateFlag = scaleFlag = rotateFlag = false;
  displayLabels = displaySamples = false;
  slicesFlag = true;
  displayFlag = true;
  bgColorWhite = bgColorChange = false;
  displayid = 0; // Display original volume by default. 

  toggleColorsFlag = false; displayNucsFlag = false; trajectoryMode = false; 
  neighborMode = false;
  //trajectoryid = 0;
  bfsMode = false;
  polygonMode = false;

  // Set initial default camera parameters.  These get re-initialized at setVol()!!!
  center[0] = 150; center[1] = 100; center[2] = 100;
  camera.towards = vec3(0,0,-1);
  camera.up = vec3(0, 1, 0);
  camera.right = vec3(1, 0, 0);
  camera.eye = center - 3 * 20 * camera.towards;
  camera.xfov = camera.yfov = 0.5;
  camera.neardist = 0.1 * 124; camera.fardist = 10 * 124;

  colormap.push_back(QColor(128,128,128));
  for(int h = 0; h < 360; h+=18) {
    QColor color;color.setHsv(h, 255, 255);
    colormap.push_back(color.toRgb());
  }

  colormap_bfs.push_back(QColor(128,128,128)); // 0
  colormap_bfs.push_back(QColor(0,255,0));   // green  1 
  colormap_bfs.push_back(QColor(255,0,255)); // magenta 2 
  colormap_bfs.push_back(QColor(0,0,255)); // blue 3
  colormap_bfs.push_back(QColor(255,0,0)); // red 4
  colormap_bfs.push_back(QColor(0,255,255)); // cyan 5
  colormap_bfs.push_back(QColor(255,255,0)); // yellow 6
}  


// Adapted from http://nehe.gamedev.net/article/replacement_for_gluperspective/21002/
// Replaces gluPerspective. Sets the frustum to perspective mode.
// fovY     - Field of vision in degrees in the y direction
// aspect   - Aspect ratio of the viewport
// zNear    - The near clipping distance
// zFar     - The far clipping distance
inline void perspectiveGL( GLdouble fovY, GLdouble aspect, GLdouble zNear, GLdouble zFar ) {
  const GLdouble pi = 3.1415926535897932384626433832795;
  GLdouble fH = tan( fovY / 360.0 * pi ) * zNear;
  GLdouble fW = fH * aspect;
  glFrustum( -fW, fW, -fH, fH, zNear, zFar );
}

void _3dView::initCameraGL() {
  // Set projection transformation
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  perspectiveGL(180.0*camera.yfov/M_PI, (GLdouble) width() /(GLdouble) height(), camera.neardist, camera.fardist);

  // Set camera transformation
  vec3 t = -(camera.towards);
  vec3& u = camera.up;
  vec3& r = camera.right;
  GLfloat  camera_matrix[16] = { r[0], 
				 u[0],  // col 1
				 t[0], 
				 0, 
				 // ------
				 r[1], 
				 u[1],  // col 2
				 t[1], 
				 0, 
				 // ------
				 r[2], 
				 u[2],  // col 3
				 t[2], 
				 0, 
				 // ------
				 0, 
				 0,     // col 4
				 0, 
				 1 };
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glMultMatrixf(camera_matrix);
  glTranslatef(-(camera.eye[0]), -(camera.eye[1]), -(camera.eye[2]));
}

// Draws slices through the image stack.
void _3dView::drawStackGL() {
  if(image_stack == NULL) return;  if(slicesFlag == false) return;

  glDisable(GL_LIGHTING);
  float w = image_stack->width(), h = image_stack->height(), d = image_stack->depth();

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glEnable(GL_TEXTURE_2D);
  // Pixelates voxel instead of smoothing. 
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST); 
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST); 

  // Display texture mapped quads that slice through th eimag evolume.
  if(image_stack->conf.use_GL_buffers)  glBindTexture(GL_TEXTURE_2D, xy_fbo->texture()); 
  else glBindTexture(GL_TEXTURE_2D, xy_texture); 
  glBegin(GL_QUADS);
  glTexCoord2f(0,1); glVertex3f(0,0,xy_zpos); glTexCoord2f(0,0); glVertex3f(0,h,xy_zpos);
  glTexCoord2f(1,0); glVertex3f(w,h,xy_zpos); glTexCoord2f(1,1); glVertex3f(w,0,xy_zpos); 
  glEnd();

  if(image_stack->conf.use_GL_buffers) glBindTexture(GL_TEXTURE_2D, yz_fbo->texture());
  else glBindTexture(GL_TEXTURE_2D, yz_texture); 
  glBegin(GL_QUADS);
  glTexCoord2f(0,1); glVertex3f(yz_xpos,0,0); glTexCoord2f(0,0); glVertex3f(yz_xpos,0,d);
  glTexCoord2f(1,0); glVertex3f(yz_xpos,h,d); glTexCoord2f(1,1); glVertex3f(yz_xpos,h,0); 
  glEnd();

  if(image_stack->conf.use_GL_buffers) glBindTexture(GL_TEXTURE_2D, xz_fbo->texture());
  else glBindTexture(GL_TEXTURE_2D, xz_texture); 
  glBegin(GL_QUADS);
  glTexCoord2f(0,1); glVertex3f(0,xz_ypos,0); glTexCoord2f(0,0); glVertex3f(0,xz_ypos,d);
  glTexCoord2f(1,0); glVertex3f(w,xz_ypos,d); glTexCoord2f(1,1); glVertex3f(w,xz_ypos,0);
  glEnd();
  glDisable(GL_TEXTURE_2D);


  // Draw lines around entire volume
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glColor3f(0,1,1); 
  glLineWidth(2);
  glBegin(GL_QUADS);
  glVertex3f(0,0,0); glVertex3f(w,0,0); glVertex3f(w,h,0); glVertex3f(0,h,0);
  glVertex3f(0,0,0); glVertex3f(0,h,0); glVertex3f(0,h,d); glVertex3f(0,0,d);
  glVertex3f(w,h,d); glVertex3f(0,h,d); glVertex3f(0,0,d); glVertex3f(w,0,d);
  glVertex3f(w,h,d); glVertex3f(w,h,0); glVertex3f(w,0,0); glVertex3f(w,0,d);
  glEnd();

  // Draw lines around orthogonal planes through the data.
  glColor3f(1,0,0); 
  glBegin(GL_QUADS);
  glVertex3f(0,0,xy_zpos); glVertex3f(0,h,xy_zpos); glVertex3f(w,h,xy_zpos); glVertex3f(w,0,xy_zpos); 
  glEnd();

  glColor3f(1,0,1); 
  glBegin(GL_QUADS);
  glVertex3f(0,xz_ypos,0); glVertex3f(0,xz_ypos,d); glVertex3f(w,xz_ypos,d); glVertex3f(w,xz_ypos,0);
  glEnd();

  glColor3f(1,1,0); 
  glBegin(GL_QUADS);
  glVertex3f(yz_xpos,0,0); glVertex3f(yz_xpos,0,d); glVertex3f(yz_xpos,h,d); glVertex3f(yz_xpos,h,0); 
  glEnd();

  // Restore color, line width, and polygon mode.
  glColor3f(1,1,1); 
  glLineWidth(1);
  glEnable(GL_LIGHTING);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}


// Draw cells that have their "displayed" flag set to true. 
void _3dView::drawCellsGL() {
  // Skip if no image stack, display is off, or edit mode is set.
  if(image_stack == NULL || editFlag || !displayFlag) return;
  GLfloat colorGL[4]; 

  if(trajectoryMode) {
    if(quickselect_traj.size() > 0) {
      /*for(size_t t = 0; t < project->stacks.size(); t += 5) {
	map<int,int> &traj_id2cell = project->stacks[t]->traj_id2cell;
	for(size_t q = 0; q < quickselect_traj.size(); q++) {
	  if(traj_id2cell.find(quickselect_traj[q]) == traj_id2cell.end()) continue;
	  // Map found, get and draw cell.
	  vector<edge::image_cell *> &cells = project->stacks[t]->cells;
	  if(q % 2 == 0) 
	    drawCellGL(cells.at(traj_id2cell[quickselect_traj[q]]), greenGL, grayGL);
	  else 
	    drawCellGL(cells.at(traj_id2cell[quickselect_traj[q]]), magentaGL, grayGL);
	}
	}*/
      GLfloat dif[4];
      dif[3] = 1;
      for(size_t q = 0; q < quickselect_traj.size(); q++) {
	map<int,int> &traj_id2cell = image_stack->traj_id2cell;
	if(traj_id2cell.find(quickselect_traj[q]) != traj_id2cell.end()) {
	  // Map found, get and draw cell.
	  vector<edge::image_cell *> &cells = image_stack->cells;
	  edge::image_cell *cur = cells.at(traj_id2cell[quickselect_traj[q]]);

	  int L = (quickselect_traj[q] % 20) + 1;
	  dif[0] = colormap[L].redF(); dif[1] = colormap[L].greenF(); dif[2] = colormap[L].blueF();
	  if(!neighborMode)  {
	    if(displayLabels && cur->mesh_face_neighbors.size() > 0)
	      drawCellLabeledNeighbors(cur, cells);
	    else 
	      drawCellGL(cur, dif, grayGL); 
	  }
	  else {
	    // Displays surface labels. 
	    if(cur->mesh_face_neighbors.size() > 0) {
	      drawCellLabeledNeighbors(cur, cells);
	    }
	    else drawCellGL(cur, dif, grayGL); 
	    vector<int> &neighbors = cur->neighbors;
	    for(size_t n = 0; n < neighbors.size(); n++) {
	      if(cells.at(neighbors[n])->trajectoryid < 0) 
		L = (cells.at(neighbors[n])->idx % 20) + 1;
	      else 
		L = (cells.at(neighbors[n])->trajectoryid % 20) + 1;
	      dif[0] = colormap[L].redF(); dif[1] = colormap[L].greenF(); dif[2] = colormap[L].blueF();
	      drawCellGL(cells.at(neighbors[n]), dif, magentaGL);
	    }
	  }

	}      
      }
    }
    /*else {
    // Get mapping from trajectory ID to cell index.
      map<int,int> &traj_id2cell = image_stack->traj_id2cell;
      if(traj_id2cell.find(trajectoryid) != traj_id2cell.end()) {
	// Map found, get and draw cell.
	vector<edge::image_cell *> &cells = image_stack->cells;
	GLfloat dif[4];
	dif[3] = 1;
	int L = (trajectoryid % 20) + 1;
	dif[0] = colormap[L].redF(); dif[1] = colormap[L].greenF(); dif[2] = colormap[L].blueF();

	if(cells.at(traj_id2cell[trajectoryid])->labels.size() > 0) {
	  drawCellLabeled(cells.at(traj_id2cell[trajectoryid]), cells);
	}
	else drawCellGL(cells.at(traj_id2cell[trajectoryid]), dif, grayGL);

	if(neighborMode) {
	  edge::image_cell *cur = cells.at(traj_id2cell[trajectoryid]);
	  vector<int> &neighbors = cur->neighbors;
	  for(size_t n = 0; n < neighbors.size(); n++) {
	    if(cells.at(neighbors[n])->trajectoryid < 0) 
	      L = (cells.at(neighbors[n])->idx % 20) + 1;
	    else 
	      L = (cells.at(neighbors[n])->trajectoryid % 20) + 1;
	    dif[0] = colormap[L].redF(); dif[1] = colormap[L].greenF(); dif[2] = colormap[L].blueF();
	    drawCellGL(cells.at(neighbors[n]), dif, magentaGL);
	  }
	}
      }
      }*/

    /*vector<edge::image_stack *> &stacks = project->stacks;
    for(size_t t = 0; t < stacks.size(); t += 4) {
      vector<edge::image_cell *> &cells = stacks[t]->cells;
      for(size_t cidx = 0; cidx < cells.size(); cidx++) {
	if(cells[cidx]->trajectoryid == trajectoryid) {
	  if(cells[cidx]->trajectoryid == trajectoryid) drawCellGL(cells[cidx], greenGL, grayGL);
	}
      }
      }*/
  }
  else if(bfsMode) {
    vector<edge::image_cell *> &cells = image_stack->cells;
    GLfloat dif[4];
    dif[3] = 1;
    for(size_t cidx = 0; cidx < cells.size(); cidx++) {
      // Draw cell only if displayed flag is set.
      if(cells[cidx]->displayed == true) {
	int L = cells[cidx]->bfs_label + 1;
	if(L > 6) L = 0;
	dif[0] = colormap_bfs[L].redF(); dif[1] = colormap_bfs[L].greenF(); dif[2] = colormap_bfs[L].blueF();
	drawCellGL(cells[cidx], dif, grayGL);
      }
    }
  }
  else {
    if(displayNucsFlag) {
      vector<edge::image_cell *> &nucs = image_stack->nuclei;
      for(int i = 0; i < 4; i++) colorGL[i] = magentaGL[i];
      for(size_t cidx = 0; cidx < nucs.size(); cidx++) {
	// Draw cell only if displayed flag is set.
	if(nucs[cidx]->displayed == true) {
	  drawCellGL(nucs[cidx], colorGL, grayGL);
	}
      }
    }
    else {
      vector<edge::image_cell *> &cells = image_stack->cells;
      for(int i = 0; i < 4; i++) colorGL[i] = greenGL[i];
      for(size_t cidx = 0; cidx < cells.size(); cidx++) {
	// Draw cell only if displayed flag is set.
	if(cells[cidx]->displayed == true) {
	  drawCellGL(cells[cidx], colorGL, grayGL);
	}
      }
    }
  }
}

// Draw selected cell with blue triangles. 
void _3dView::drawCellGL(edge::image_cell *cell, GLfloat front[4], GLfloat back[4], bool skip_gl_buffers) {
  //vector<vec3> &v0 = cell->v0, &v1 = cell->v1;

  /*  glDisable(GL_LIGHTING);  
  glLineWidth(3);
  glBegin(GL_LINES);
  // Draw any boundaries present in the cell as lines.
  for(size_t i = 0; i < v0.size(); i++) {
    glVertex3f(v0[i][0], v0[i][1], v0[i][2]); glVertex3f(v1[i][0], v1[i][1], v1[i][2]);
  }
  glEnd();
  glLineWidth(1);
  glEnable(GL_LIGHTING);*/

  // Draw primary axes.
  /*glDisable(GL_LIGHTING);  
  glLineWidth(3);
  glBegin(GL_LINES);
  vec3 &c = cell->mesh->centroid;
  for(size_t i = 0; i < 3; i++) {
    glVertex3f(c[0], c[1], c[2]);
    vec3 x = c + 1.2 * cell->mesh->axes3[i];
    glVertex3f(x[0], x[1], x[2]);
  }
  glEnd();
  glLineWidth(1); // Restore line width.*/


  /*glDisable(GL_LIGHTING);  
  glLineWidth(3);
  glBegin(GL_LINES);
  vec3 &c = cell->mesh->centroid;
  for(size_t i = 1; i < 3; i++) {
    glVertex3f(c[0], c[1], c[2]);
    vec3 x = c + 1.2 * cell->mesh->axes3[i];
    glVertex3f(x[0], x[1], x[2]);
  }
  glEnd();
  glLineWidth(1);*/ // Restore line width.


  /*{
  glLineWidth(3);
  glBegin(GL_LINES);
  vec3 &c = cell->cp1;
  for(size_t i = 0; i < 3; i++) {
    glVertex3f(c[0], c[1], c[2]);
    vec3 x = c + 1.2 * cell->mesh->axes3[i];
    glVertex3f(x[0], x[1], x[2]);
  }
  glEnd();
  glLineWidth(1); // Restore line width.
  }
  {
  glLineWidth(3);
  glBegin(GL_LINES);
  vec3 &c = cell->cp2;
  for(size_t i = 0; i < 3; i++) {
    glVertex3f(c[0], c[1], c[2]);
    vec3 x = c + 1.2 * cell->mesh->axes3[i];
    glVertex3f(x[0], x[1], x[2]);
  }
  glEnd();
  glLineWidth(1); // Restore line width.
  }*/

  //  glEnable(GL_LIGHTING);
  if(polygonMode) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glLineWidth(2);
  }

  // Draw vertex buffer data for mesh.  TODO: Handle triangle colors as well.
  glMaterialfv(GL_BACK, GL_DIFFUSE,  back); glMaterialfv(GL_FRONT, GL_DIFFUSE, front);

  if(image_stack->conf.use_GL_buffers && skip_gl_buffers == false) {
    //glBindBufferARB(GL_ARRAY_BUFFER_ARB, cell->mesh->vboId);
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, cell->mesh->vboId);

    glEnableClientState(GL_NORMAL_ARRAY);
    //glEnableClientState(GL_COLOR_ARRAY);
    glEnableClientState(GL_VERTEX_ARRAY);

    glNormalPointer(GL_FLOAT, 0, (void*)(cell->num_faces() * 3 * sizeof(vec3)));
    //glColorPointer(3, GL_FLOAT, sizeof(triangle), (void*)(sizeof(vec3) * 4));
    glVertexPointer(3, GL_FLOAT, 0, 0);

    glDrawArrays(GL_TRIANGLES, 0, 3*cell->num_faces());

    glDisableClientState(GL_VERTEX_ARRAY); 
    //glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);

    glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
  }
  else {
    vector<triangle> &tribuf = cell->mesh->tribuf;    
    glBegin(GL_TRIANGLES);
    for(size_t i = 0; i < tribuf.size(); i++) {
      triangle &t = tribuf[i];
      glNormal3f(t.normal[0], t.normal[1], t.normal[2]);
      glVertex3f(t[2][0], t[2][1], t[2][2]);
      glVertex3f(t[1][0], t[1][1], t[1][2]);
      glVertex3f(t[0][0], t[0][1], t[0][2]);
    }
    glEnd();
    
  }

  if(polygonMode) {
    glLineWidth(2);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  }

}

void _3dView::drawSelectedGL() {
  if(image_stack == NULL || selected_cell == NULL || !displayFlag) return;
  if(editFlag && displayLabels && selected_cell->labels.size() > 0) {
    // Draw cell part or draw labels.
    if(selected_cell->parts.size() > 0) {
      vector<edge::image_cell *> &parts = selected_cell->parts;
      for(size_t p = 0; p < parts.size(); p++) {
	if(parts[p]->displayed) {
	  drawCellGL(parts[p], magentaGL, grayGL, true);
	}
      }
      if(selected_part != NULL) drawCellGL(selected_part, redGL, grayGL, true);
    }
    else drawCellLabeled(selected_cell);
  }
  else {
    if(selected_cell->v_samples.size() > 0) drawCellLabeled(selected_cell);
    else drawCellGL(selected_cell, blueGL, grayGL);
  }
}

// TODO: Push labeling into mesh object. Use color portion of VBOs. 
void _3dView::drawCellLabeled(edge::image_cell *cell) {
  /*vector<vec3> &v0 = cell->v0, &v1 = cell->v1;
  glDisable(GL_LIGHTING);  
  glLineWidth(3);
  glBegin(GL_LINES);
  // Draw any boundaries present in the cell as lines.
  for(size_t i = 0; i < v0.size(); i++) {
    glVertex3f(v0[i][0], v0[i][1], v0[i][2]); glVertex3f(v1[i][0], v1[i][1], v1[i][2]);
  }
  glEnd();
  glLineWidth(1);*/

  glEnable(GL_LIGHTING);
  if(!displaySamples) {
    // Displays surface labels. 
    if(cell->labels.size() > 0) {
      GLfloat dif[4];
      vector<int> &labels = cell->labels[cell->label_idx];    
      vector<triangle> &tribuf = cell->mesh->tribuf;
      dif[3] = 1;
      for(size_t i = 0; i < labels.size(); i++) {
	int L = labels[i] == 0 ? 0 : (((labels[i] - 1) % 20) + 1);
	dif[0] = colormap[L].redF(); dif[1] = colormap[L].greenF(); dif[2] = colormap[L].blueF();
	glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
	glBegin(GL_TRIANGLES);
	triangle &t = tribuf[i];
	glNormal3f(t.normal[0], t.normal[1], t.normal[2]);
	glVertex3f(t[2][0], t[2][1], t[2][2]);
	glVertex3f(t[1][0], t[1][1], t[1][2]);
	glVertex3f(t[0][0], t[0][1], t[0][2]);
	glEnd();
      }
    }
  }
  else {
    // Displays surface projection.
    GLfloat dif[4];
    vector<triangle> &tribuf = cell->mesh->tribuf;
    vector<uint8> &samples = cell->v_samples;
    for(size_t i = 0; i < samples.size(); i++) {
      float L = samples[i];
      float val = L/255.0f;
      dif[0] = val; dif[1] = val; dif[2] = 0; dif[3] = 1;
      glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
      glBegin(GL_TRIANGLES);
      triangle &t = tribuf[i];
      glNormal3f(t.normal[0], t.normal[1], t.normal[2]);
      glVertex3f(t[2][0], t[2][1], t[2][2]);
      glVertex3f(t[1][0], t[1][1], t[1][2]);
      glVertex3f(t[0][0], t[0][1], t[0][2]);
      glEnd();
    }
  }
}

void _3dView::drawCellLabeledNeighbors(edge::image_cell *cell, vector<edge::image_cell *> &cells) {
  // Displays surface labels. 
  if(cell->mesh_face_neighbors.size() > 0) {
    GLfloat dif[4];
    vector<int> &labels = cell->mesh_face_neighbors;    
    vector<triangle> &tribuf = cell->mesh->tribuf;
    dif[3] = 1;
    for(size_t i = 0; i < labels.size(); i++) {
      if(labels[i] < 0) {
	dif[0] = colormap[0].redF(); dif[1] = colormap[0].greenF(); dif[2] = colormap[0].blueF();
      }
      else {
	int traj_id = cells.at(labels[i])->trajectoryid;
	if(traj_id < 0) traj_id = cells.at(labels[i])->idx;
	int L = (traj_id % 20) + 1;
	dif[0] = colormap[L].redF(); dif[1] = colormap[L].greenF(); dif[2] = colormap[L].blueF();
      }
      glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
      glBegin(GL_TRIANGLES);
      triangle &t = tribuf[i];
      glNormal3f(t.normal[0], t.normal[1], t.normal[2]);
      glVertex3f(t[2][0], t[2][1], t[2][2]);
      glVertex3f(t[1][0], t[1][1], t[1][2]);
      glVertex3f(t[0][0], t[0][1], t[0][2]);
      glEnd();
    }
  }
}


void _3dView::initializeGL() {
  // Initialize graphics modes 
  glShadeModel(GL_SMOOTH);
  glEnable(GL_MULTISAMPLE); 
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_NORMALIZE);  
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  glEnable(GL_SMOOTH); 
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);  // Black clear color, initially.
  // NOTE: Back-face rendering is disabled!!!
  // glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
}

void _3dView::paintGL() {
  if(bgColorChange) {
    if(bgColorWhite) glClearColor(1.0f, 1.0f, 1.0f, 0.0f); else glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    bgColorChange = false;
  }
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  initCameraGL();   // Initialize the OpenGL camera model.
  drawStackGL();    // Draw orthogonal slices from image volume/stack.

  if(image_stack != NULL) {
    drawSelectedGL(); // Draw selected mesh/cell.
    drawCellsGL();    // Draw non-selected meshes/cells.
  }

  // Create a "head lamp" as the only light source.
  glLoadIdentity();
  GLfloat light_position[] = { 0.0, 0.0, -1.0, 0.0 };
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glEnable(GL_LIGHT0);  
}

void _3dView::resizeGL(int width, int height) {
  // Resize window
  glViewport(0, 0, width, height);
  // Resize camera vertical field of view to match aspect ratio of viewport
  camera.yfov = atanf(tanf(camera.xfov) * (float) height/ (float) width); 
  updateGL();
}

void _3dView::mousePressEvent(QMouseEvent *event) { 
  if(image_stack == NULL) return;
  lastPos = event->pos(); lastPosFlag = true; 

  if(undisplayFlag) {
    // Create a ray from position where mouse was released.  Use this
    // ray to turn off the display of a cell.
    vec3 o, p; //camera.create_ray(o, p, event->x(), height() - event->y() - 1, width(), height());
    camera.create_ray(o, p, event->x(), event->y(), width(), height());
    edge::image_cell *cell = NULL;
    if(editFlag && selected_cell != NULL) cell = selected_cell->select_part(o, p); else cell = image_stack->select_cell(o, p);
    if(cell != NULL) cell->displayed = false;
    undisplayFlag = false;
    event->accept(); updateGL();
  }
  // TODO: State information is screwed up especially with multiple time points.
  if(selectFlag) {
    // Use this ray to a "select" an individual cell. 
    vec3 o, p;
    camera.create_ray(o, p, event->x(), event->y(), width(), height());
    if(editFlag && selected_cell != NULL) {
      edge::image_cell *cell = selected_cell->select_part(o, p);
      if(cell != NULL) {    
	// TODO: Have some way to unselect a selected cell.
	if(selected_part != NULL) selected_part->displayed = true; 
	selected_part = cell;
	selected_part->displayed = false;
	updateMeasurementTree(selected_part);
      }
      else if(selected_part != NULL) { selected_part->displayed = true; selected_part = NULL; }
    }
    else {
      edge::image_cell *cell = image_stack->select_cell(o, p);
      if(cell != NULL) { 
	// TODO: Have some way to unselect a selected cell.
	if(selected_cell != NULL) {
	  selected_cell->displayed = true; 
	  selected_cell->selected = false;
	}
	selected_cell = cell;

	// TODO: Allow selection of multiple cells. 
	// if(event->modifier() == Qt::ControlModifier) { }

	// This new cell is "undisplayed" since it is handled separately by 3dview. 
	// Not sure if this is the best decision.
	selected_cell->displayed = false;
	selected_cell->selected = true;
	updateMeasurementTree(selected_cell);
      } // Nothing we selected, unselect cell. 
      else if(selected_cell != NULL) { 
	updateMeasurementTree(selected_cell);
	selected_cell->selected = false;
	selected_cell->displayed = true; 
	selected_cell = NULL; 
      }
    }
    selectFlag = false;
    event->accept(); updateGL();
  }
}

void _3dView::mouseMoveEvent(QMouseEvent *event) {
  if(image_stack == NULL) return;
  // Record a last position if none has been set.  
  if(!lastPosFlag) { lastPos = event->pos(); lastPosFlag = true; return; }
  int dx = event->x() - lastPos.x(), dy = -(event->y() - lastPos.y());
  lastPos = event->pos(); // Remember mouse position.
  if(dx == 0 &&  dy == 0) return; // No movement case.

  // Capture paint event before changing the scene.
  /*if(paintFlag && image_stack != NULL && selected_cell != NULL && (event->buttons() & Qt::LeftButton)) {
    vec3 o, p;
    camera.create_ray(o, p, event->x(), event->y(), width(), height());
    selected_cell->label(o, p, paintLabel);
    event->accept(); updateGL();
  }
  // Otherwise, change the scene position by rotation , scaling, or translation.
  else*/ 
  if (rotateFlag || (event->buttons() & Qt::LeftButton)) { // Rotate scene around a center point.
    camera.rotate_scene(dx, dy, width(), height(), center); event->accept(); updateGL();
  } 
  else if (scaleFlag || (event->buttons() & Qt::MidButton)) { // Scale the scene.
    camera.scale_scene(dx, dy, width(), height(), center); event->accept(); updateGL();
  }
  else if (translateFlag || (event->buttons() & Qt::RightButton)) { // Translate the scene.
    camera.translate_scene(dx, dy, width(), height(), center); event->accept(); updateGL();
  }
}

void _3dView::addItem2Tree(QString key, QString value, QString units) {
  QTreeWidgetItem *w = new QTreeWidgetItem(QStringList() << key << value << units);
  zTree->addTopLevelItem(w);
}

void _3dView::updateItem(int itemNum, QString value) { zTree->topLevelItem(itemNum)->setText(1, value); }

// Updates measurements as cells are clicked.
void _3dView::updateMeasurementTree(edge::image_cell *cell) {
  if(zTree == NULL || cell == NULL) return;
  QTreeWidgetItem *w = NULL;
  int analysis_id = cell->analysis_id();
  if(analysis_id == analysis::CellShapes ||
     analysis_id == analysis::NucsOnly ||
     analysis_id == analysis::NucsMembranes ) {
    if(zTree->topLevelItemCount() == 0) {
      addItem2Tree("Traj. ID", QString::number(cell->trajectoryid), "um2");
      addItem2Tree("Triangles", QString::number(cell->num_faces()), "");
      addItem2Tree("Volume", QString::number(cell->model_volume), "um3");
      addItem2Tree("Surface Area", QString::number(cell->surface_area), "um2");
      addItem2Tree("Cell Index", QString::number(cell->idx + 1), "");
      addItem2Tree("Bend", QString::number(cell->max_cell_bend), "au");
      addItem2Tree("Nuc. Apical Dist.", QString::number(cell->nuc_apical_dist), "um");
      vector<float> &pca_dims = cell->pca_dims;
      float anisotropy = pca_dims.at(1) / pca_dims.at(0);
      addItem2Tree("anisotropy", QString::number(anisotropy), "");      
    }
    else {
      updateItem(0, QString::number(cell->trajectoryid));
      updateItem(1, QString::number(cell->num_faces()));
      updateItem(2, QString::number(cell->model_volume));
      updateItem(3, QString::number(cell->surface_area));
      updateItem(4, QString::number(cell->idx + 1));
      updateItem(5, QString::number(cell->max_cell_bend));
      updateItem(6, QString::number(cell->nuc_apical_dist));
      vector<float> &pca_dims = cell->pca_dims;
      float anisotropy = pca_dims.at(1) / pca_dims.at(0);
      updateItem(7, QString::number(anisotropy));      
    }
  }
  else if(cell->analysis_id() == analysis::DorsalFolds) {
    edge::dorsalfolds_t &M = cell->dorsalfolds;
    if(zTree->topLevelItemCount() == 0) {
      // TODO: Replace these with function calls.
      addItem2Tree("Traj. ID", QString::number(cell->trajectoryid),"");
      addItem2Tree("Cell Index", QString::number(cell->idx + 1), "");
      addItem2Tree("Triangles", QString::number(cell->num_faces()), "");
      w = new QTreeWidgetItem(QStringList() << tr("Total Surf Area") << QString::number(M.total_surface_area) << tr("um2"));
      zTree->addTopLevelItem(w);
      w = new QTreeWidgetItem(QStringList() << tr("Vol") << QString::number(M.total_volume) << tr("um3"));
      zTree->addTopLevelItem(w);
      w = new QTreeWidgetItem(QStringList() << tr("Above Baz Vol") << QString::number(M.above_baz_volume) << tr("um3"));
      zTree->addTopLevelItem(w);
      w = new QTreeWidgetItem(QStringList() << tr("Above Baz Surf Area") << QString::number(M.above_baz_surface_area) << tr("um2"));
      zTree->addTopLevelItem(w);
      w = new QTreeWidgetItem(QStringList() << tr("Centroid X") << QString::number(M.cent_x) << tr(""));
      zTree->addTopLevelItem(w);
      w = new QTreeWidgetItem(QStringList() << tr("Centroid Y") << QString::number(M.cent_y) << tr(""));
      zTree->addTopLevelItem(w);
      w = new QTreeWidgetItem(QStringList() << tr("Centroid Z") << QString::number(M.cent_z) << tr(""));
      zTree->addTopLevelItem(w);
      w = new QTreeWidgetItem(QStringList() << tr("Total Length") << QString::number(M.total_length) << tr("um"));
      zTree->addTopLevelItem(w);
      w = new QTreeWidgetItem(QStringList() << tr("Apical Length") << QString::number(M.apical_length) << tr("um"));
      zTree->addTopLevelItem(w);
      addItem2Tree("Baz Position", QString::number(M.bazooka_pos), "%");
      addItem2Tree("Total Par Intensity", QString::number(M.total_par_intensity), "");
      addItem2Tree("Basal Par Intensity", QString::number(M.basal_par_intensity), "");
      addItem2Tree("Baz Patch Vol", QString::number(M.bazooka_patch_volume),"um3");
      addItem2Tree("Baz Patch Intensity", QString::number(M.bazooka_patch_intensity),"");
      addItem2Tree("Baz Non-patch Vol", QString::number(M.bazooka_nonpatch_volume),"um3");
      addItem2Tree("Baz Non-patch Intensity", QString::number(M.bazooka_nonpatch_intensity),"");
      addItem2Tree("Basal Par 4-voxel Vol", QString::number(M.basal_4voxel_par_volume),"um3");
      addItem2Tree("Basal Par 4-voxel Intensity", QString::number(M.basal_4voxel_par_intensity),"");
      addItem2Tree("Basal Par 2-voxel Vol", QString::number(M.basal_2voxel_par_volume),"um3");
      addItem2Tree("Basal Par 2-voxel Intensity", QString::number(M.basal_2voxel_par_intensity),"");
    }
    else {
      updateItem(0, QString::number(cell->trajectoryid));
      updateItem(1, QString::number(cell->idx + 1));
      updateItem(2, QString::number(cell->num_faces()));
      zTree->topLevelItem(3)->setText(1, QString::number(M.total_surface_area));
      zTree->topLevelItem(4)->setText(1, QString::number(M.total_volume));
      zTree->topLevelItem(5)->setText(1, QString::number(M.above_baz_volume));
      zTree->topLevelItem(6)->setText(1, QString::number(M.above_baz_surface_area));
      zTree->topLevelItem(7)->setText(1, QString::number(M.cent_x));
      zTree->topLevelItem(8)->setText(1, QString::number(M.cent_y));
      zTree->topLevelItem(9)->setText(1, QString::number(M.cent_z));
      zTree->topLevelItem(10)->setText(1, QString::number(M.total_length));
      zTree->topLevelItem(11)->setText(1, QString::number(M.apical_length));
      zTree->topLevelItem(12)->setText(1, QString::number(M.bazooka_pos));
      zTree->topLevelItem(13)->setText(1, QString::number(M.total_par_intensity));
      zTree->topLevelItem(14)->setText(1, QString::number(M.basal_par_intensity));
      zTree->topLevelItem(15)->setText(1, QString::number(M.bazooka_patch_volume));
      zTree->topLevelItem(16)->setText(1, QString::number(M.bazooka_patch_intensity));
      zTree->topLevelItem(17)->setText(1, QString::number(M.bazooka_nonpatch_volume));
      zTree->topLevelItem(18)->setText(1, QString::number(M.bazooka_nonpatch_intensity));
      zTree->topLevelItem(19)->setText(1, QString::number(M.basal_4voxel_par_volume));
      zTree->topLevelItem(20)->setText(1, QString::number(M.basal_4voxel_par_intensity));
      zTree->topLevelItem(21)->setText(1, QString::number(M.basal_2voxel_par_volume));
      zTree->topLevelItem(22)->setText(1, QString::number(M.basal_2voxel_par_intensity));
    }
  }
  zTree->resizeColumnToContents(0);
}

