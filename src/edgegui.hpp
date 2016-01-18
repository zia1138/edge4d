#ifndef __EDGEGUI_HPP__
#define __EDGEGUI_HPP__

#define GL_GLEXT_PROTOTYPES
#include <QtGui>
#include <QHash>

#include "ui_edge.h"
#include "ui_param.h"
#include "ui_quickselect.h"
#include "ui_measure.h"

#include "edge.hpp"

class EdgeWin : public QMainWindow, private Ui::EdgeWin {
  Q_OBJECT
public:
  EdgeWin(QMainWindow *parent = 0);
  ~EdgeWin() { if(project != NULL) delete project; }
protected:
  // TODO: Create a GUI parameter set and a separate processing parameter parameter set!
  edge::hyperstack *project;
  int cur_stack;
  QHash<QTreeWidgetItem *, int> mTree2Idx; 
  void keyPressEvent(QKeyEvent *event);
  void progressDialog(QString message, vector<QThread *> threads);

  QString getFileName(QString postfix);

  // Updates info display below sliders.
  void updateInfo() {
    if(project != NULL) {
      edgeconfig &conf = project->conf;
      QString info;
      info += "Frame: " + QString::number(cur_stack * conf.stepframe  + conf.startframe + 1) + " ";
      info += "Time: " + QString::number((cur_stack * conf.stepframe + conf.startframe) * conf.time_step) + " sec. ";
      if(view != NULL) {
	double alpha, x,y,z; 
	alpha = project->conf.voxel_alpha();
	view->getPos(x,y,z);
	info += "Slider: (" + QString::number(x) + ", " + QString::number(y) + ", " + QString::number(z) + ") voxels; ";
	info += "(" + QString::number(alpha*x) + ", " + QString::number(alpha*y) + ", " + QString::number(alpha*z) + ") microns";
      }
      infoLabel->setText(info);
    }
  }
  void load_process_manual_groundtruth(string filename);
  void updateVolumes() {
    int vx = scrollXplane->value(), vy = scrollYplane->value(), vz = scrollZplane->value();
    view->setYZplane(vx); view->setXZplane(vy); view->setXYplane(vz);
  }

private slots:
  // Menu triggers
  void on_actionLoad_Stack_triggered(); // Load data
  void on_action_Quit_triggered() { if(project != NULL) delete project; exit(0); } // Quit.
  void on_actionSave_Labels_triggered(); // Save .raw segmentation file for MATLAB analysis.
  void on_actionSave_Ortho_Slices_triggered(); // Saves .PNGs slicing through data.
  void on_actionSave_Processing_triggered(); // Saves .PNGS of display volumes.

  void on_actionSave_Analysis_triggered();

  void on_actionSave_Trajectory_Lengths_triggered();
  void on_actionSave_Stable_Neighbors_triggered();
  void on_actionSave_Neighbor_Swaps_triggered();
  void on_actionSave_Neighbor_Analysis_triggered();
  void on_actionSave_Neighbor_Counts_triggered();
  void on_actionSave_All_Trajectories_triggered();
  void on_actionSelect_1_Neighbors_triggered();
  void on_actionSelect_2_Neighbors_triggered();
  void on_actionSelect_1_2_Neighbors_Right_triggered();
  void on_actionSelect_1_2_Neighbors_Left_triggered();
  void on_actionSelect_Custom_Neighbors_All_triggered();

  void on_actionSave_Contact_Analysis_triggered();

  void on_actionToggle_Colors_triggered() { if(view != NULL) { view->toggleColors(); view->updateGL(); } }
  void on_actionToggle_Nuclei_triggered() { 
    if(project == NULL || view == NULL) return;
    view->toggleDisplayNucs(); 
    if(project->conf.analysis_id == analysis::ACME_Compare || 
       project->conf.analysis_id == analysis::StaticAnalysis) {
      project->stacks[cur_stack]->toggle_all_nucs();
    }
    else {
      project->stacks[cur_stack]->match_nuclei_to_cells_displayed();
    }
    view->updateGL(); 
  }

  //  Breaks meshes up into parts (submeshes).
  void on_action_Segment_triggered() {
    if(project == NULL) return;
    edge::image_cell *selected = view->getSelectedCell();
    if(selected != NULL && selected->parts.size() == 0) {
      selected->split_parts();
      QTreeWidgetItem *model = mTree->topLevelItem(selected->idx);
      vector<edge::image_cell *> &parts = selected->parts;
      for(size_t i = 0; i < parts.size(); i++) {
	QTreeWidgetItem *part = new QTreeWidgetItem(model, QStringList() << "Part " + QString::number(parts[i]->idx + 1));
	part->setSelected(true);
      }
    }
  }

  void on_action_Measure_triggered();

  void on_actionUpdate_Measurement_Parameters_triggered();

  void on_updateButton_clicked(bool); 
  void on_refreshButton_clicked(bool); 
  void on_actionSelect_Cell_Shapes_triggered();
  void on_actionSelect_Trajectories_triggered();

  void on_actionMeasurements_triggered() { if(!zDock->isVisible()) zDock->show();  }
  void on_actionModels_triggered() { if(!mDock->isVisible()) mDock->show(); }

  // Scroll (x,y,z,t).
  void on_scrollXplane_valueChanged(int); void on_scrollYplane_valueChanged(int);
  void on_scrollZplane_valueChanged(int); void on_scrollTime_valueChanged(int);
};


// Parameter dialog box.
class ParamWin : public QDialog, private Ui::ParamWin {
  Q_OBJECT
public:
  ParamWin(edgeconfig &conf_, QWidget *parent = 0);
protected:
  edgeconfig &conf;
  double voxel2microns(int voxels);   // Convert voxel distance to microns.
  double voxel2microns3(int voxels);  // Convert voxel volume to microns^3. 
  void updateMicrons();
  void updateBoundaryValues(int idx);
private slots:
  void on_loadButton_clicked(bool);
  void on_cancelButton_clicked(bool) { reject(); }  // User hit cancel.
  // Update micron conversion.
  void on_spinVoxelX_valueChanged(double v) { spinVoxelY->setValue(v); updateMicrons(); }
  void on_spinVoxelZ_valueChanged(double)  { updateMicrons(); }
  void on_checkScaleToZ_stateChanged(int)  { updateMicrons(); }
  void on_spinScale_valueChanged(double)   { updateMicrons(); }
  void on_spinR1_valueChanged(int)         { updateMicrons(); }
  void on_spinMaxHoleRad_valueChanged(int) { updateMicrons(); }
  void on_spinMaxCompVol_valueChanged(int) { updateMicrons(); }
  void on_spinMinCompVol_valueChanged(int) { updateMicrons(); }
  void on_spinErrorVoxThresh_valueChanged(int) { updateMicrons(); }
  void on_listBoundaryRanges_currentRowChanged(int row) { updateBoundaryValues(row); }
};

// Quick select a subset of cells.
class QuickSelectDialog : public QDialog, private Ui::QuickSelectDialog {
  Q_OBJECT
public:
  QuickSelectDialog(QWidget *parent = 0) : QDialog(parent) { setupUi(this); }  
  void getSelected(vector<int> &selected) {
    selected.clear();
    QStringList list = textSelect->toPlainText().split(QRegExp("(\\s+|,)"));
    for(int i = 0; i < list.size(); i++) {
      bool ok;
      int number = list[i].toInt(&ok);
      if(ok) selected.push_back(number);
    }
  }
  void setText(QString txt) {  textSelect->setText(txt); }
private slots:
};

// Displays and allows adjustment measurement configuration parameters.
class MeasureWin : public QDialog, private Ui::MeasureWin {
  Q_OBJECT
public:
  MeasureWin(edgeconfig &conf, QWidget *parent = 0);
protected:
  edgeconfig &conf;
private slots:  
  void accept();
};


#endif
