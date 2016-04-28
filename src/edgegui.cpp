#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

#include <math.h>

#include "edgegui.hpp"

using namespace std;

EdgeWin::EdgeWin(QMainWindow *parent) : QMainWindow(parent) { 
  setupUi(this); 
  project = NULL; 
  cur_stack = 0;
  view->setMeasurementTree(zTree);
}

// Displays a "busy wait" progress dialog
void EdgeWin::progressDialog(QString message, vector<QThread *> threads) {
  // Start all of the threads.
  for(int i = 0; i < (int)threads.size(); i++) threads[i]->start();
  /*QProgressDialog progress(message, 0, 0, 0, this); // if start = end = 0, then busy wait
  progress.setWindowTitle("EDGE4D Data Processing");
  progress.setModal(true); 
  progress.show();*/
  //QApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
  while(true) {
    // Wait until both loader threads complete.
    bool done = true;
    for(size_t i = 0; i < threads.size(); i++) { 
      done = done && threads[i]->wait(100);// QApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    if(done) break;
  }
}

// Updates 3dview to match selcted items in cell list.
void EdgeWin::on_updateButton_clicked(bool) {
  if(project == NULL) return;
  //cout << "update displayed cells" << endl;
  for(int i = 0; i < project->num_cells(cur_stack); i++) {
    project->stacks[cur_stack]->cells[i]->displayed = mTree->topLevelItem(i)->isSelected();
    QTreeWidgetItem *model = mTree->topLevelItem(i);
    for(int k = 0; k < model->childCount(); k++) {
      project->stacks[cur_stack]->cells[i]->parts[k]->displayed = model->child(k)->isSelected();
    }
  }
  view->updateGL();
}
// Make cell list match cells displayed in 3d view.
void EdgeWin::on_refreshButton_clicked(bool) {
  // TODO: Use selectedItems in QTreeWidget!!!
  if(project == NULL) return;
  for(int i = 0; i < project->num_cells(cur_stack); i++) {
    QTreeWidgetItem *model = mTree->topLevelItem(i);
    model->setSelected(project->stacks[cur_stack]->cells[i]->displayed);
    for(int k = 0; k < model->childCount(); k++) {
      model->child(k)->setSelected(project->stacks[cur_stack]->cells[i]->parts[k]->displayed);
    }
  }
}

void EdgeWin::on_actionSelect_Cell_Shapes_triggered() {
  if(project == NULL) return;
  QuickSelectDialog quickselect(this);
  if(quickselect.exec() == QDialog::Rejected) return; // User hit cancel
  vector<int> selected;
  quickselect.getSelected(selected);
  int ncells = project->num_cells(cur_stack);
  // Undisplay all cells.
  for(int i = 0; i < ncells; i++) {
    project->stacks[cur_stack]->cells[i]->displayed = false;
  }
  for(size_t s = 0; s < selected.size(); s++) {
    if(selected[s] > 0 && selected[s] <= ncells) {
      project->stacks[cur_stack]->cells.at(selected[s]-1)->displayed = true;
    }
  }
  view->updateGL();
}



void EdgeWin::on_actionSelect_Trajectories_triggered() {
  if(project == NULL || view == NULL) return;
  QuickSelectDialog quickselect(this);

  vector<int> cur_selected_traj;
  view->getQuickSelectTraj(cur_selected_traj);
  sort(cur_selected_traj.begin(), cur_selected_traj.end());

  QString selected_traj_str;
  if(cur_selected_traj.size() > 0) {
    selected_traj_str = QString::number(cur_selected_traj[0]);
    for(size_t c = 1; c < cur_selected_traj.size(); c++) 
      selected_traj_str += ", " + QString::number(cur_selected_traj[c]);
  }
  quickselect.setText(selected_traj_str);

  if(quickselect.exec() == QDialog::Rejected) return; // User hit cancel
  vector<int> selected;
  quickselect.getSelected(selected);
  view->setQuickSelectTraj(selected);
  view->updateGL();
}

void EdgeWin::on_actionSelect_1_Neighbors_triggered() {
  if(project == NULL || view == NULL) return;
  vector<int> quick_sel, neighbor_sel;
  view->getQuickSelectTraj(quick_sel);
  project->stacks[cur_stack]->get_neighbors_traj(neighbor_sel, 1, quick_sel);
  view->setQuickSelectTraj(neighbor_sel);
  view->updateGL();
}

void EdgeWin::on_actionSelect_2_Neighbors_triggered() {
  if(project == NULL || view == NULL) return;
  vector<int> quick_sel, neighbor_sel;
  view->getQuickSelectTraj(quick_sel);
  project->stacks[cur_stack]->get_neighbors_traj(neighbor_sel, 2, quick_sel);
  view->setQuickSelectTraj(neighbor_sel);
  view->updateGL();
}

void EdgeWin::on_actionSelect_1_2_Neighbors_Right_triggered() {
  if(project == NULL || view == NULL) return;
  vector<int> quick_sel, neighbor_sel;
  view->getQuickSelectTraj(quick_sel);
  vector<int> order;
  QuickSelectDialog quickselect(this);
  if(quickselect.exec() == QDialog::Rejected) return; // User hit cancel
  quickselect.getSelected(order);
  if(order.size() >  0) {
    project->stacks[cur_stack]->get_neighbors_traj(neighbor_sel, edge::Neighbor_Xgreater, order, quick_sel);
    view->setQuickSelectTraj(neighbor_sel);
  }
  view->updateGL();
}

void EdgeWin::on_actionSelect_Custom_Neighbors_All_triggered(){
  if(project == NULL || view == NULL) return;
  vector<int> quick_sel, neighbor_sel;
  view->getQuickSelectTraj(quick_sel);
  vector<int> order;
  QuickSelectDialog quickselect(this);
  if(quickselect.exec() == QDialog::Rejected) return; // User hit cancel
  quickselect.getSelected(order);
  if(order.size() >  0) {
    project->stacks[cur_stack]->get_neighbors_traj(neighbor_sel, edge::Neighbor_All, order, quick_sel);
    view->setQuickSelectTraj(neighbor_sel);
    view->updateGL();
  }
}

void EdgeWin::on_actionSelect_1_2_Neighbors_Left_triggered() {
  if(project == NULL || view == NULL) return;
  vector<int> quick_sel, neighbor_sel;
  view->getQuickSelectTraj(quick_sel);
  vector<int> order;
  QuickSelectDialog quickselect(this);
  if(quickselect.exec() == QDialog::Rejected) return; // User hit cancel
  quickselect.getSelected(order);
  if(order.size() > 0) {
    project->stacks[cur_stack]->get_neighbors_traj(neighbor_sel, edge::Neighbor_Xless, order, quick_sel);
    view->setQuickSelectTraj(neighbor_sel);
  }
  view->updateGL();
}

// ShortcutA keyboard interface.
void EdgeWin::keyPressEvent(QKeyEvent *event ) {
  if(view != NULL) {
    switch(event->key()) {
      // Cycle through display volumes. 
    case Qt::Key_D: view->toggleDOG();       break;
    case Qt::Key_C: view->toggleDOGrev();    break;
    case Qt::Key_F: view->toggleDisplay();   break;
      // Set rotate, scale, translate modes (for track pad mainly).
    case Qt::Key_R: view->toggleRotate();    break;
    case Qt::Key_S: view->toggleScale();     break;
    case Qt::Key_T: view->toggleTranslate(); break;
      // Select and undislay by selection.
    case Qt::Key_X: view->toggleUndisplay(); break;
    case Qt::Key_A: view->toggleSelected();  break;
    case Qt::Key_E: view->toggleEdit();      break;
    case Qt::Key_M: view->toggleSlices();    break;
    case Qt::Key_W: view->toggleBgColor();   break;
      // Surface projection and segmentation of models
    case Qt::Key_L: view->toggleLabels();   break;
    case Qt::Key_7: view->cycleLabels();    break;
    case Qt::Key_8: view->toggleSamples();  break;
    case Qt::Key_G: view->toggleTrajectoryMode(); break;
    case Qt::Key_N: view->toggleNeighborMode(); break;
    case Qt::Key_B: view->toggleBfsMode(); break;
    case Qt::Key_P: view->togglePolygonMode(); break;
    case Qt::Key_Plus:
    case Qt::Key_Equal: view->add_selected_to_quickselect(); break;
    case Qt::Key_Minus: view->remove_selected_from_quickselect(); break;
    default: break;
    }
    event->accept();
    view->updateGL();
  }
  else QMainWindow::keyPressEvent(event);
}

void EdgeWin::on_actionSave_Analysis_triggered() {
  if(project == NULL) return;
  QString defaultname = project->basename.c_str();
  if(project->conf.analysis_id == analysis::StaticAnalysis )  
    defaultname += "_static_analysis.txt";
  else if(project->conf.analysis_id == analysis::ManualGroundTruth)
    defaultname += "_manual_groundtruth.txt";
  else if(project->conf.analysis_id == analysis::DorsalFolds) 
    defaultname += "_dorsalfolds.txt";
  else
    return;

  QString fileName = QFileDialog::getSaveFileName(this, tr("Save Tab Delimiated"), defaultname, tr("TXT (*.txt)"));
  if(fileName == "") return;
  if(project->conf.analysis_id == analysis::ManualGroundTruth) 
    load_process_manual_groundtruth(fileName.toStdString());
  else if(project->conf.analysis_id == analysis::DorsalFolds) {
    vector<int> cur_selected_traj;
    view->getQuickSelectTraj(cur_selected_traj);
    if(cur_selected_traj.size() == 0) 
      project->save_dorsalfolds(fileName.toStdString());
    else
      project->save_dorsalfolds(fileName.toStdString(), cur_selected_traj);
  }
  else
    project->global_analysis(fileName.toStdString());
  
}

// Gets filename for a given analaysis.
QString EdgeWin::getFileName(QString postfix) {
  if(project == NULL) return QString("");
  QString defaultname = project->basename.c_str();
  defaultname += postfix;
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save Tab Delimiated"), defaultname, tr("TXT (*.txt)"));
  return fileName;
}

void EdgeWin::on_actionSave_Trajectory_Lengths_triggered() {
  QString fileName = getFileName("_trajectory_lengths.txt"); if(fileName == "") return;
  project->save_trajectory_lengths(fileName.toStdString());
}

void EdgeWin::on_actionSave_Stable_Neighbors_triggered() {
  QString fileName = getFileName("_stable_neighbors.txt"); if(fileName == "") return;
  project->save_stable_neighbors(fileName.toStdString());
}

void EdgeWin::on_actionSave_Neighbor_Swaps_triggered() {
  QString fileName = getFileName("_neighbor_swaps.txt"); if(fileName == "") return;
  project->save_neighbor_swaps(fileName.toStdString());
}

void EdgeWin::on_actionSave_Neighbor_Analysis_triggered() {
  QString fileName = getFileName("_neighbor_analysis.txt"); if(fileName == "") return;
  project->save_neighbor_analysis(fileName.toStdString());
}

void EdgeWin::on_actionSave_Neighbor_Counts_triggered() {
  QString fileName = getFileName("_neighbor_counts.txt"); if(fileName == "") return;
  project->save_neighbor_counts(fileName.toStdString());
}

struct MeasureThread : public QThread {
  edge::hyperstack *project;
  MeasureThread(edge::hyperstack *project) : project(project) { }
  void run() { project->measure_trajectories_threaded();  }
};

void EdgeWin::on_action_Measure_triggered() {
  if(project == NULL) return;
  MeasureThread measure_thread(project);
  vector<QThread *> threads; threads.push_back((QThread *)&measure_thread);
  progressDialog("Measuring. Please wait.", threads);
  view->updateMeasurementTree();
  updateVolumes();
  view->updateGL();
}

struct SaveAllTrajThread : public QThread {
  int analysis_id;
  edge::hyperstack *project;
  QString fileName;
  SaveAllTrajThread(int analysis_id, edge::hyperstack *project, QString fileName) : 
    analysis_id(analysis_id), project(project), fileName(fileName) { }
  void run() { 
    switch(analysis_id) {
    case 0: project->save_analyze_trajectories(fileName.toStdString()); break;
    case 1: project->save_contact_analysis(fileName.toStdString()); break;
    default: break;
    }
  }
};

void EdgeWin::on_actionSave_All_Trajectories_triggered() {
  QString fileName = getFileName("_all_traj.txt"); if(fileName == "") return;
  vector<QThread *> threads;  
  SaveAllTrajThread save_thread(0, project, fileName);
  threads.push_back((QThread*)&save_thread);
  progressDialog("Measuring and saving trajectories. Please wait.", threads);
}


void EdgeWin::on_actionSave_Contact_Analysis_triggered() {
  QString fileName = getFileName("_contact_analysis.txt"); if(fileName == "") return;
  vector<QThread *> threads;  
  SaveAllTrajThread save_thread(1, project, fileName);
  threads.push_back((QThread*)&save_thread);
  progressDialog("Measuring and saving trajectories. Please wait.", threads);
}


// Saves a .raw image file that can be loaded in matlab to view an unprocessed segmentation.
void EdgeWin::on_actionSave_Labels_triggered() {
  if(project == NULL) return;
  QString defaultname = project->basename.c_str();
  defaultname += ".raw";
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save RAW segementation data"), defaultname, tr("RAW (*.raw)"));
  if(fileName == "") return;
  edge::image_stack *stack = project->stacks[cur_stack];

  stack->save_raw_segmentation(fileName.toStdString(), true);

  string matlabFile = 
    QFileInfo(fileName).path().toStdString() + "/" + 
    QFileInfo(fileName).baseName().toStdString() + ".m";
  
  QString matlab_commands;
  QTextStream(&matlab_commands) << 
    "% isotropic voxel dimensions" << endl <<
    "alpha = " << stack->conf.voxel_alpha() << "; % microns" << endl << 
    "w = " << stack->iwidth << ";" << endl <<
    "h = " << stack->iheight << ";" << endl << 
    "d = " << stack->idepth << ";" << endl << 
    "fid = fopen('" << fileName << "');" << endl <<
    "label_matrix = fread(fid, w * h * d, '*int32');" << endl <<
    "labels3d = reshape(label_matrix, [w h d]);" << endl <<
    "fclose(fid);" << endl;

  ofstream matlab(matlabFile.c_str());
  matlab << matlab_commands.toStdString();
  matlab.close();

  QMessageBox msgBox;
  msgBox.setIcon(QMessageBox::Information);
  msgBox.setWindowTitle("Additional Information");
  msgBox.setText("Segmentation data has been saved in raw binary format.");
  msgBox.setInformativeText("MATLAB commands for loading have been saved to " + QString(matlabFile.c_str()));

  msgBox.setIcon(QMessageBox::Information);
  msgBox.setDetailedText(matlab_commands);
  msgBox.exec();
}


void EdgeWin::on_actionSave_Processing_triggered() {
  if(project == NULL) return;

  QString basename = project->basename.c_str();
  vector<volume8 *> &dispvols = project->stacks[cur_stack]->dispvols;

  for(size_t d = 0; d < dispvols.size(); d++) {
    QString info_str;
    QTextStream info(&info_str);

    volume8 &vol = *dispvols[d];
    QImage xy_img(vol.width, vol.height, QImage::Format_RGB32);
    QImage xz_img(vol.width, vol.depth, QImage::Format_RGB32);    
    QImage yz_img(vol.height, vol.depth, QImage::Format_RGB32);    

    int vx = scrollXplane->value(), vy = scrollYplane->value(), vz = scrollZplane->value();
    //int xpos = vol.width / 2, ypos = vol.height / 2, zpos = vol.depth / 2;
    int xpos = vx, ypos = vy, zpos = vz;

    for(int x = 0; x < vol.width; x++) 
      for(int y = 0;  y < vol.height; y++) { uint8 v = vol(x,y,zpos); xy_img.setPixel(x,y, qRgb(v,v,v)); }

    for(int x = 0; x < vol.width; x++) 
      for(int z = 0; z < vol.depth; z++) { uint8 v = vol(x,ypos,z); xz_img.setPixel(x,z, qRgb(v,v,v)); }
    
    for(int y = 0; y < vol.height; y++) 
      for(int z = 0; z < vol.depth; z++) { uint8 v = vol(xpos,y,z); yz_img.setPixel(y,z, qRgb(v,v,v)); }

    QString basename_c;
    if(d < 10) basename_c = basename + "_dispvol_0" + QString::number(d);
    else basename_c = basename + "_dispvol_" + QString::number(d);

    QString xy_file = basename_c + "_xy.png", yz_file = basename_c + "_yz.png", xz_file = basename_c + "_xz.png";

    info << xy_file << endl << yz_file << endl << xz_file << endl;
    xy_img.save(xy_file, "PNG"); yz_img.save(yz_file, "PNG"); xz_img.save(xz_file, "PNG");

  }
  
}

void EdgeWin::on_actionSave_Ortho_Slices_triggered() {
  if(project == NULL) return;

  QString basename = project->basename.c_str();

  QString info_str;
  QTextStream info(&info_str);

  volume8 *channels[2];
  channels[0] = project->stacks[cur_stack]->edge_channel();
  channels[1] = project->stacks[cur_stack]->nuc_channel();

  basename += "_f" + QString::number(project->conf.startframe + cur_stack);

  // Build images with orthogonal slice data.
  for(int c = 0; c < 2; c++) {
    if(channels[c] == NULL) continue;

    volume8 &vol = *channels[c];
    QImage xy_img(vol.width, vol.height, QImage::Format_RGB32);
    QImage xz_img(vol.width, vol.depth, QImage::Format_RGB32);    
    QImage yz_img(vol.height, vol.depth, QImage::Format_RGB32);    

    int vx = scrollXplane->value(), vy = scrollYplane->value(), vz = scrollZplane->value();
    //int xpos = vol.width / 2, ypos = vol.height / 2, zpos = vol.depth / 2;
    int xpos = vx, ypos = vy, zpos = vz;

    for(int x = 0; x < vol.width; x++) 
      for(int y = 0;  y < vol.height; y++) { uint8 v = vol(x,y,zpos); xy_img.setPixel(x,y, qRgb(v,v,v)); }

    for(int x = 0; x < vol.width; x++) 
      for(int z = 0; z < vol.depth; z++) { uint8 v = vol(x,ypos,z); xz_img.setPixel(x,z, qRgb(v,v,v)); }
    
    for(int y = 0; y < vol.height; y++) 
      for(int z = 0; z < vol.depth; z++) { uint8 v = vol(xpos,y,z); yz_img.setPixel(y,z, qRgb(v,v,v)); }

    QString basename_c = basename + "_c" + QString::number(c);
    QString xy_file = basename_c + "_xy.png", yz_file = basename_c + "_yz.png", xz_file = basename_c + "_xz.png";

    info << xy_file << endl << yz_file << endl << xz_file << endl;
    xy_img.save(xy_file, "PNG"); yz_img.save(yz_file, "PNG"); xz_img.save(xz_file, "PNG");

  }

  // Build a 2-color image to save.
  if(channels[0] != NULL && channels[1] != NULL) {
    volume8 &vol1 = *channels[0];
    volume8 &vol2 = *channels[1];

    QImage xy_img(vol1.width, vol1.height, QImage::Format_RGB32);
    QImage xz_img(vol1.width, vol1.depth, QImage::Format_RGB32);    
    QImage yz_img(vol1.height, vol1.depth, QImage::Format_RGB32);    

    int xpos = vol1.width / 2, ypos = vol1.height / 2, zpos = vol1.depth / 2;

    for(int x = 0; x < vol1.width; x++) 
      for(int y = 0;  y < vol1.height; y++) { 
	uint8 v1 = vol1(x,y,zpos); 
	uint8 v2 = vol2(x,y,zpos); 
	xy_img.setPixel(x,y, qRgb(v2,v1,v2)); 
      }

    for(int x = 0; x < vol1.width; x++) 
      for(int z = 0; z < vol1.depth; z++) { 
	uint8 v1 = vol1(x,ypos,z); 
	uint8 v2 = vol2(x,ypos,z); 
	xz_img.setPixel(x,z, qRgb(v2,v1,v2)); 
      }
    
    for(int y = 0; y < vol1.height; y++) 
      for(int z = 0; z < vol1.depth; z++) { 
	uint8 v1 = vol1(xpos,y,z); 
	uint8 v2 = vol2(xpos,y,z); 
	yz_img.setPixel(y,z, qRgb(v2,v1,v2)); 
      }

    QString basename_c = basename + "_combined" ;
    QString xy_file = basename_c + "_xy.png", yz_file = basename_c + "_yz.png", xz_file = basename_c + "_xz.png";

    info << xy_file << endl << yz_file << endl << xz_file << endl;
    xy_img.save(xy_file, "PNG"); yz_img.save(yz_file, "PNG"); xz_img.save(xz_file, "PNG");
  }


  QMessageBox msgBox;
  msgBox.setIcon(QMessageBox::Information);
  msgBox.setWindowTitle("Additional Information");
  msgBox.setText("Orthogonal slices have been saved to .PNG files.");
  msgBox.setDetailedText(info_str);
  msgBox.exec();
}

// Compares a 2-d segmentation region obtained from orthogonally slicing the data
// to a ground truth region generated by hand image annotation.
float compare_groundtruth(float voxel_alpha, vector<ivec3> &truth, vector<ivec3> &test) {
  vector<ivec3> combined;
  for(size_t t = 0; t < truth.size(); t++) combined.push_back(truth[t]);
  for(size_t s = 0; s < test.size(); s++) combined.push_back(test[s]);

  // Get bounding box around combined voxel set.
  ivec3 mn, mx;
  geom::bounding_box(mn, mx, combined);
  ivec3 p0(mn[0], mn[1], 0);
  int w = (mx[0] - mn[0] + 1) + 2, h = (mx[1] - mn[1] + 1) + 2;

  // Shift so voxels fit in bounding region.
  for(size_t i = 0; i < truth.size(); i++) truth[i] = truth[i] - p0 + ivec3(1,1,0);
  for(size_t i = 0; i < test.size(); i++) test[i] = test[i] - p0 + ivec3(1,1,0);    

  // Get outlines of the ground truth and test shapes
  vector<ivec3> truth_outline, test_outline;
  volume8 outline(w,h);

  outline.fill(0); 
  outline.set(truth, 255); 
  outline.fill_holes();
  outline.dilate(1, 255, 0);
  outline.dilate(1, 0, 255);
  outline.outline_vox(truth_outline);

  outline.fill(0); 
  outline.set(test, 255); 
  outline.fill_holes();
  outline.dilate(1, 255, 0);
  outline.dilate(1, 0, 255);
  outline.outline_vox(test_outline);

  // Compute EDT of the truth outline data.
  outline.fill(0);
  for(size_t i = 0; i < truth_outline.size(); i++) outline(truth_outline[i]) = 255;
  volume32 EDTsq(outline.width, outline.height, outline.depth);
  EDTsq.ComputeEDT(outline);

  // Use EDT to comptue RMS between the shapes in microns.
  float N = test_outline.size();
  float RMS = 0;
  for(size_t i = 0; i < test_outline.size(); i++) {
    // Convert distance to microns using voxel_alpha. 
    float dist_microns = float(voxel_alpha) * sqrtf(EDTsq(test_outline[i]));
    RMS += dist_microns * dist_microns / N;
  }
  RMS = sqrtf(RMS);

  // Restore position of test and truth voxels.
  for(size_t i = 0; i < truth.size(); i++) truth[i] = truth[i] + p0 - ivec3(1,1,0);  
  for(size_t i = 0; i < test.size(); i++) test[i] = test[i] + p0 - ivec3(1,1,0);;

  return RMS;
}

// Loads image files with ground truth data.
void EdgeWin::load_process_manual_groundtruth(string filename) {
  if(project == NULL) return;
  QString basename = project->basename.c_str();   // Scan for ground truth files.

  basename += "_f" + QString::number(project->conf.startframe + cur_stack);

  volume8 *channels[2] = { NULL, NULL };
  channels[0] =  project->stacks[cur_stack]->edge_channel();
  channels[1] = project->stacks[cur_stack]->nuc_channel();

  // Check to see if all ground truth files are present.
  QString info_str;
  QTextStream info(&info_str);
  int missing = 0;
  for(int c = 0; c < 2; c++) {
    if(channels[c] == NULL) continue;
    QString basename_c = basename + "_c" + QString::number(c);
    // NOTE: They must be appended with " copy.png".
    QString xy_file = basename_c + "_xy copy.png";
    QString yz_file = basename_c + "_yz copy.png";
    QString xz_file = basename_c + "_xz copy.png";
    if(!QFile::exists(xy_file)) { info << xy_file << endl; missing++; }
    if(!QFile::exists(yz_file)) { info << yz_file << endl; missing++; }
    if(!QFile::exists(xz_file)) { info << xz_file << endl; missing++; }
  }
  if(missing > 0) {
    QMessageBox msgBox;
    msgBox.setIcon(QMessageBox::Information);
    msgBox.setWindowTitle("Files Missing");
    msgBox.setText("Unable to find the following ground truth files.");
    msgBox.setDetailedText(info_str);
    msgBox.exec();    
    return;
  }

  ofstream output(filename.c_str());
  output << "data.id\t";  output << "seg.size\t"; output << "status\t"; 
  output << "RMS" << endl;

  // All ground truth files found. Now evaluate.
  for(int c = 0; c < 2; c++) {
    if(channels[c] == NULL) continue;
    volume8 &vol = *channels[c];

    QString basename_c = basename + "_c" + QString::number(c);
    QString xy_file = basename_c + "_xy copy.png";
    QString yz_file = basename_c + "_yz copy.png";
    QString xz_file = basename_c + "_xz copy.png";

    QImage xy_img(xy_file, "PNG"), xz_img(xz_file, "PNG"), yz_img(yz_file, "PNG");    

    // Check dimensions of ground truth files.
    if(xy_img.width() != vol.width  || xy_img.height() != vol.height ||
       yz_img.width() != vol.height || yz_img.height() != vol.depth ||
       xz_img.width() != vol.width  || xz_img.height() != vol.depth) { 
      QMessageBox msgBox;
      msgBox.setWindowTitle("Ground Truth Processing Error");
      msgBox.setText("Invalid dimension found in ground truth files.");
      msgBox.exec();
      return;
    }

    // Get ground truth slices.
    volume8 xy_gt(vol.width, vol.height), yz_gt(vol.height, vol.depth), xz_gt(vol.width, vol.depth);

    for(int x = 0; x < xy_img.width(); x++) 
      for(int y = 0; y < xy_img.height(); y++) {  
	QColor v(xy_img.pixel(x,y));
	if(v.red() == 0 && v.green() == 255 && v.blue() == 0) xy_gt(x,y) = 255;
      }

    for(int y = 0; y < yz_img.width(); y++) 
      for(int z = 0; z < yz_img.height(); z++) {  
	QColor v(yz_img.pixel(y,z));
	if(v.red() == 0 && v.green() == 255 && v.blue() == 0) yz_gt(y,z) = 255;
      }

    for(int x = 0; x < xz_img.width(); x++) 
      for(int z = 0; z < xz_img.height(); z++) {  
	QColor v(xz_img.pixel(x,z));
	if(v.red() == 0 && v.green() == 255 && v.blue() == 0) xz_gt(x,z) = 255;
      }

    // Get ground truth regions.
    vector< vector<ivec3> > comps_xy, comps_yz, comps_xz;
    xy_gt.components(comps_xy, 0, 0); yz_gt.components(comps_yz, 0, 0); xz_gt.components(comps_xz, 0, 0);

    // Clear largest. It's the background.
    geom::ClearLargest(comps_xy); geom::ClearLargest(comps_yz); geom::ClearLargest(comps_xz);

    // Generate segmentation images and compare.
    int xpos = vol.width / 2, ypos = vol.height / 2, zpos = vol.depth / 2;

    // Slice through each region and obtain segmented components
    vector< vector<ivec3> > xy_seg, yz_seg, xz_seg;
    if(c == 0) { // For cell segmetnation.
      project->stacks[cur_stack]->ortho_comps_cells(yz_seg, 0, xpos);
      project->stacks[cur_stack]->ortho_comps_cells(xz_seg, 1, ypos);
      project->stacks[cur_stack]->ortho_comps_cells(xy_seg, 2, zpos);
    }
    else if(c == 1) {  // For nuclei.
      project->stacks[cur_stack]->ortho_comps_nucs(yz_seg, 0, xpos);
      project->stacks[cur_stack]->ortho_comps_nucs(xz_seg, 1, ypos);
      project->stacks[cur_stack]->ortho_comps_nucs(xy_seg, 2, zpos);
    }
    else exit(1); 

    //  For each slice, transfer 2d portion to an (x,y,0) image.
    for(size_t s = 0; s < xy_seg.size(); s++) {
      vector<ivec3> &seg = xy_seg[s];
      for(size_t i = 0; i < seg.size(); i++) seg[i].z() = 0;
    }
    // Transfer (x,z) dimensions to (x,y,0).
    for(size_t s = 0; s < xz_seg.size(); s++) {
      vector<ivec3> &seg = xz_seg[s];
      for(size_t i = 0; i < seg.size(); i++) { seg[i].y() = seg[i].z();seg[i].z() = 0; }
    }
    // Transfer (y,z) to (x,y,0)
    for(size_t s = 0; s < yz_seg.size(); s++) {
      vector<ivec3> &seg = yz_seg[s];
      for(size_t i = 0; i < seg.size(); i++) { seg[i].x() = seg[i].y(); seg[i].y() = seg[i].z(); seg[i].z() = 0;  }
    }

    // TODO: Figure out if 0.2 overlap % should be a parameter.
    // Also clean up repeated code. 

    int missed = 0;
    vector<int> seg_status, assigned_truth;
    edge::static_analysis_internal(seg_status, assigned_truth, 
				   missed, ivec3(vol.width, vol.height, 1),
				   0.2, 0.2, comps_xy, xy_seg);

    output << "c" << c << "_xy\t" << missed << "\t-1\t-2" << endl;// channel and dim

    for(size_t a = 0; a < assigned_truth.size(); a++) {
      if(xy_seg[a].size() == 0) continue;
      output << "c" << c << "_xy" << '\t'; // channel and dim
      output << xy_seg[a].size() << '\t'; // 
      output << seg_status[a] << '\t';  // status
      if(assigned_truth[a] < 0) output << "-1" << endl;
      else {
	vector<ivec3> &seg_reg =   xy_seg[a];
	vector<ivec3> &truth_reg = comps_xy.at(assigned_truth[a]);
	// Compare these regions.
	float RMS = compare_groundtruth(project->conf.voxel_alpha(), seg_reg, truth_reg);
	output << RMS << endl;
      }
    }

    
    edge::static_analysis_internal(seg_status, assigned_truth, 
				   missed, ivec3(vol.width, vol.depth, 1),
				   0.2, 0.2, comps_xz, xz_seg);

    output << "c" << c << "_xz\t" << missed << "\t-1\t-2" << endl;// channel and dim

    for(size_t a = 0; a < assigned_truth.size(); a++) {
      if(xz_seg[a].size() == 0 ) continue;
      output << "c" << c << "_xz" << '\t'; // channel and dim
      output << xz_seg[a].size() << '\t';
      output << seg_status[a] << '\t';
      if(assigned_truth[a] < 0) output << "-1" << endl;
      else {
	vector<ivec3> &seg_reg =   xz_seg[a];
	vector<ivec3> &truth_reg = comps_xz.at(assigned_truth[a]);
	// Compare these regions.
	float RMS = compare_groundtruth(project->conf.voxel_alpha(), seg_reg, truth_reg);
	output << RMS << endl;
      }
    }

    edge::static_analysis_internal(seg_status, assigned_truth, 
				   missed, ivec3(vol.height, vol.depth, 1),
				   0.2, 0.2, comps_yz, yz_seg);

    output << "c" << c << "_yz\t" << missed << "\t-1\t-2" << endl;// channel and dim

    for(size_t a = 0; a < assigned_truth.size(); a++) {
      if(yz_seg[a].size() == 0) continue;
      output << "c" << c << "_yz" << '\t'; // channel and dim
      output << yz_seg[a].size() << '\t';
      output << seg_status[a] << '\t';
      if(assigned_truth[a] < 0) output << "-1" << endl;
      else {
	vector<ivec3> &seg_reg =   yz_seg[a];
	vector<ivec3> &truth_reg = comps_yz.at(assigned_truth[a]);
	// Compare these regions.
	float RMS = compare_groundtruth(project->conf.voxel_alpha(), seg_reg, truth_reg);
	output << RMS << endl;
      }
    }
    // Advance to next channel.
  }
  
}


// Adjust orthogonal slices when scroll bars are changed.
void EdgeWin::on_scrollXplane_valueChanged(int value) { if(view != NULL) { view->setYZplane(value); view->updateGL(); updateInfo(); } }
void EdgeWin::on_scrollYplane_valueChanged(int value) { if(view != NULL) { view->setXZplane(value); view->updateGL(); updateInfo(); } }
void EdgeWin::on_scrollZplane_valueChanged(int value) { if(view != NULL) { view->setXYplane(value); view->updateGL(); updateInfo(); } }

// Change displayed volume when slider is advanced.
void EdgeWin::on_scrollTime_valueChanged(int value) {
  //cout << "time(value)=" << value << endl;
  if(view != NULL) {
    cur_stack = value; 
    //cout << "update measurement info when time bar is changed!!" << endl;
    mTree->clear(); 
    mTree2Idx.clear();
    // TODO: Organize into noise, background, etc.
    // at time point 0 get models 
    for(int idx = 0; idx < project->num_cells(cur_stack); idx++) {
      QString defaultText = tr("Model ") + QString::number(idx + 1);
      QTreeWidgetItem *modelItem = new QTreeWidgetItem(mTree, QStringList() << defaultText);
      modelItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
      mTree2Idx[modelItem] = idx;
    }
    view->updateVol(project->stacks[cur_stack]);
    view->updateGL();
    // Update plane slider bars to stay within current stack's size.
    int width  = project->stacks[cur_stack]->width(), 
      height = project->stacks[cur_stack]->height(), 
      depth = project->stacks[cur_stack]->depth(); 
    scrollXplane->setRange(0, width - 1);  
    scrollYplane->setRange(0, height - 1); 
    scrollZplane->setRange(0, depth - 1);  

    updateInfo();
  }
}

double ParamWin::voxel2microns(int voxels) {
  double scalef   = spinScale->value();
  double voxel_x = spinVoxelX->value() / scalef;
  double voxel_y = spinVoxelY->value() / scalef;
  double voxel_z = spinVoxelZ->value() / scalef;
  double voxel_xy = min(voxel_x, voxel_y);
  if(checkScaleToZ->checkState() == Qt::Checked) return voxel_z * voxels;
  else return voxel_xy * voxels;
}

double ParamWin::voxel2microns3(int voxels) {
  double scalef = spinScale->value();
  double voxel_x = spinVoxelX->value() / scalef;
  double voxel_y = spinVoxelY->value() / scalef;
  double voxel_z = spinVoxelZ->value() / scalef;
  double voxel_xy = min(voxel_x, voxel_y);
  if(checkScaleToZ->checkState() == Qt::Checked) return (voxel_z * voxel_z * voxel_z) * voxels;
  else return (voxel_xy * voxel_xy * voxel_xy)  * voxels;
}

void ParamWin::updateMicrons() {
  spinR1Microns->setValue(voxel2microns(spinR1->value()));
  spinMaxHoleRadMicrons->setValue(voxel2microns(spinMaxHoleRad->value()));
  spinMinCompVolMicrons->setValue(voxel2microns3(spinMinCompVol->value()));
  spinMaxCompVolMicrons->setValue(voxel2microns3(spinMaxCompVol->value()));
}

void ParamWin::updateBoundaryValues(int idx) {
  if(0 <= idx && idx < (int)conf.boundaryconf.size()) {
    boundaryconf_t &bconf = conf.boundaryconf.at(idx);
    spinBoundaryStart->setValue(bconf.startframe);
    spinBoundaryEnd->setValue(bconf.endframe);
    spinThreshHigh->setValue(bconf.thresh_high); 
    spinSigmaHigh->setValue(bconf.sigma_high);
    spinRRepair->setValue(bconf.r_repair);  
    spinR2Repair->setValue(bconf.r2_repair);
    spinBoundaryFmin->setValue(bconf.boundary_fmin);
    spinBoundaryMinComp->setValue(bconf.boundary_mincomp);
  }
}

// Initialize paramater dialog with current parameter set.
ParamWin::ParamWin(edgeconfig &conf_, QWidget *parent) : 
  QDialog(parent), conf(conf_) {
  setupUi(this);

  // system parameters
  comboAnalysis->setCurrentIndex(conf.analysis_id);
  spinThreads->setValue(conf.threads);
  spinHSThreads->setValue(conf.hs_threads);
  if(conf.run_refinement) checkRefine->setCheckState(Qt::Checked); else checkRefine->setCheckState(Qt::Unchecked);
  if(conf.run_mesh) checkRunMesh->setCheckState(Qt::Checked); else checkRunMesh->setCheckState(Qt::Unchecked);
  if(conf.keep_processing) checkKeepProcessing->setCheckState(Qt::Checked); else checkKeepProcessing->setCheckState(Qt::Unchecked);
  if(conf.keep_edge_channel_only) checkKeepEdgeOnly->setCheckState(Qt::Checked); 
  else checkKeepEdgeOnly->setCheckState(Qt::Unchecked);
  if(conf.use_GL_buffers) checkUseGLBuffers->setCheckState(Qt::Checked);
  else checkUseGLBuffers->setCheckState(Qt::Unchecked);

  // TODO: Need to tie maximum value of edge channel to spinChannels!

  spinEdgeChannel->setMaximum(conf.channels); spinEdgeChannel->setValue(conf.edge_channel);
  spinNucChannel->setMaximum(conf.channels); spinNucChannel->setValue(conf.nuc_channel);

  // tiff parameters
  spinFrames->setValue(conf.frames);
  spinSlices->setValue(conf.slices);
  spinChannels->setValue(conf.channels);
  if(conf.override_tiff) checkOverrideTIFF->setCheckState(Qt::Checked); 
  else checkOverrideTIFF->setCheckState(Qt::Unchecked);
  spinStartFrame->setValue(conf.startframe+1);
  spinEndFrame->setValue(conf.endframe+1);
  spinStepFrame->setValue(conf.stepframe);
  spinVoxelX->setValue(conf.voxelX);
  spinVoxelY->setValue(conf.voxelY);
  spinVoxelZ->setValue(conf.voxelZ);
  spinScale->setValue(conf.scalef);

  spinStartFrame->setMinimum(1); spinEndFrame->setMaximum(conf.frames);
  spinEndFrame->setMinimum(1); spinEndFrame->setMaximum(conf.frames);

  if(conf.scaleToZ) checkScaleToZ->setCheckState(Qt::Checked); 
  else checkScaleToZ->setCheckState(Qt::Unchecked);

  spinTimeStep->setValue(conf.time_step);

  // enhance parameters
  spinHistNx->setValue(conf.hist_nx);
  spinHistNy->setValue(conf.hist_ny);
  spinHistNz->setValue(conf.hist_nz);
  spinFmin->setValue(conf.fmin);
  spinFmax->setValue(conf.fmax);

  spinR1->setValue(conf.r1); spinR2->setValue(conf.r2);
  spinAlpha1->setValue(conf.alpha1); spinAlpha2->setValue(conf.alpha2);

  spinR1_Nuc->setValue(conf.r1_nuc); spinR2_Nuc->setValue(conf.r2_nuc);
  spinAlpha1_Nuc->setValue(conf.alpha1_nuc); spinAlpha2_Nuc->setValue(conf.alpha2_nuc);

  spinRNonEdge->setValue(conf.r_nonedge);

  checkRepairNucs->setCheckState(conf.repair_nucs ? Qt::Checked : Qt::Unchecked);
  
  updateBoundaryValues(0);

  // DOG parameters
  for(size_t b = 1; b < conf.boundaryconf.size(); b++) {
    QString rangeStr = "[" + QString::number(conf.boundaryconf[b].startframe) + ", " + 
      QString::number(conf.boundaryconf[b].endframe) + "]";
    listBoundaryRanges->addItem(rangeStr);
  }

  spinThreshLow->setValue(conf.thresh_low); spinThreshNuc->setValue(conf.thresh_nuc);
  spinSigmaLow->setValue(conf.sigma_low);   spinSigmaNuc->setValue(conf.sigma_nuc);

  spinMinNucVol->setValue(conf.min_nuc_vol); spinMaxNucVol->setValue(conf.max_nuc_vol);


  // segmentation parameters
  spinMaxHoleRad->setValue(conf.max_hole_rad);
  spinMaxCompVol->setValue(conf.max_comp_vol); 
  spinMinCompVol->setValue(conf.min_comp_vol);
  spinNoiseCompVol->setValue(conf.noise_comp_vol);
  spinInternalLimit->setValue(conf.internal_limit);

  spinMCdim->setValue(conf.mcdim);
  spinRefineSigma->setValue(conf.refine_sigma);
  spinRefineDistance->setValue(conf.refine_dist);
  spinRefineAlpha->setValue(conf.refine_alpha);
  spinRefineStep->setValue(conf.refine_stepsize);
  spinRefineIter->setValue(conf.refine_iter);

  spinNeighborContactAlpha->setValue(conf.neighbor_alpha);
}

// Collect user updated parameter values from dialog entries. 
void ParamWin::on_loadButton_clicked(bool) { 
  // system parameters
  conf.analysis_id = comboAnalysis->currentIndex();
  conf.threads = spinThreads->value();
  conf.hs_threads = spinHSThreads->value();
  conf.run_refinement = checkRefine->checkState() == Qt::Checked;
  conf.run_mesh = checkRunMesh->checkState() == Qt::Checked;
  conf.keep_processing = checkKeepProcessing->checkState() == Qt::Checked;
  conf.keep_edge_channel_only = checkKeepEdgeOnly->checkState() == Qt::Checked;
  conf.edge_channel = spinEdgeChannel->value();
  conf.nuc_channel = spinNucChannel->value();
  conf.use_GL_buffers = checkUseGLBuffers->checkState() == Qt::Checked;
  
  // tiff parameters
  conf.frames = spinFrames->value();
  conf.slices = spinSlices->value();
  conf.channels = spinChannels->value();
  conf.override_tiff = checkOverrideTIFF->checkState() == Qt::Checked;
  conf.startframe = spinStartFrame->value() - 1;
  conf.stepframe = spinStepFrame->value();
  conf.endframe = spinEndFrame->value() - 1;

  if(conf.startframe > conf.endframe) conf.endframe = conf.startframe;

  conf.time_step = spinTimeStep->value();

  conf.voxelX = spinVoxelX->value();
  conf.voxelY = spinVoxelY->value();
  // NOTE: This is critical for backward compatibility.
  conf.voxelXY = min(conf.voxelX, conf.voxelY);
  conf.voxelZ = spinVoxelZ->value();
  conf.scalef = spinScale->value();
  conf.scaleToZ = checkScaleToZ->checkState() == Qt::Checked;
  
  // enhance parameters
  conf.hist_nx = spinHistNx->value();
  conf.hist_ny = spinHistNy->value();
  conf.hist_nz = spinHistNz->value();
  conf.fmin = spinFmin->value();
  conf.fmax = spinFmax->value();
  conf.r1 = spinR1->value(); conf.r2 = spinR2->value();
  conf.alpha1 = spinAlpha1->value(); conf.alpha2 = spinAlpha2->value();

  conf.r1_nuc = spinR1_Nuc->value(); conf.r2_nuc = spinR2_Nuc->value();
  conf.alpha1_nuc = spinAlpha1_Nuc->value(); conf.alpha2_nuc = spinAlpha2_Nuc->value();

  conf.r_nonedge = spinRNonEdge->value();
  conf.repair_nucs = checkRepairNucs->checkState() == Qt::Checked;

  // dog parameters
  int row = listBoundaryRanges->currentRow();
  if(row < 0) row = 0; // udpate default
  boundaryconf_t &bconf = conf.boundaryconf.at(row);
  bconf.boundary_fmin = spinBoundaryFmin->value();
  bconf.boundary_mincomp = spinBoundaryMinComp->value();
  bconf.sigma_high = spinSigmaHigh->value();  
  bconf.thresh_high = spinThreshHigh->value(); 
  bconf.r_repair = spinRRepair->value();
  bconf.r2_repair = spinR2Repair->value();

  conf.sigma_low = spinSigmaLow->value();   conf.sigma_nuc = spinSigmaNuc->value();
  conf.thresh_low = spinThreshLow->value(); conf.thresh_nuc = spinThreshNuc->value();

  // mesh parameters
  conf.min_nuc_vol = spinMinNucVol->value();
  conf.max_nuc_vol = spinMaxNucVol->value();

  conf.max_hole_rad = spinMaxHoleRad->value();
  conf.max_comp_vol = spinMaxCompVol->value(); 
  conf.min_comp_vol = spinMinCompVol->value();
  conf.internal_limit = spinInternalLimit->value();
  conf.noise_comp_vol = spinNoiseCompVol->value();

  conf.neighbor_alpha = spinNeighborContactAlpha->value();

  conf.mcdim = spinMCdim->value();
  conf.refine_sigma = spinRefineSigma->value();
  conf.refine_sigma = spinRefineSigma->value();
  conf.refine_dist = spinRefineDistance->value();
  conf.refine_alpha = spinRefineAlpha->value();
  conf.refine_stepsize = spinRefineStep->value();
  conf.refine_iter = spinRefineIter->value();

  accept(); 
}

struct ProcessThread : public QThread {
  edge::hyperstack *project;
  ProcessThread(edge::hyperstack *project) : project(project) { }
  void run() { project->process(); }
};

void EdgeWin::on_actionLoad_Stack_triggered() {
  QString stackFile = QFileDialog::getOpenFileName(this, "Open File", "", "TIFF (*.tif *.tiff)");
  if(stackFile == "") return;
  edgeconfig conf;

  // Load XML configuration file if it exists.
  QString suffix = QFileInfo(stackFile).completeSuffix();
  QString baseName = QFileInfo(stackFile).baseName();
  QString xmlfile = stackFile;
  xmlfile.replace("." + suffix, ""); xmlfile += ".xml";
  if(QFile(xmlfile).exists()) conf.load_xml(xmlfile.toStdString());

  // Load parameters from TIFF.
  if(conf.override_tiff == false) conf.load_tiff(stackFile.toStdString());

  // Now, display configuration dialog. 
  ParamWin configure(conf, this); 
  if(configure.exec() == QDialog::Rejected) return; // User hit cancel.

  QTime time; time.start();
  
  conf.save_xml(xmlfile.toStdString());   // Save XML file with configuration information.

  edge::hyperstack *old_project = NULL;
  if(project != NULL) {  old_project = project; cur_stack = 0;  }

  project = new edge::hyperstack(stackFile.toStdString(), conf);
  project->basename = QFileInfo(stackFile).baseName().toStdString();

  vector<QThread *> threads;  ProcessThread process_thread(project);
  threads.push_back((QThread*)&process_thread);
  progressDialog("Processing image stacks. Please wait.", threads);

  if(project->num_stacks() == 0) {
    delete project;
    project = NULL;
    QMessageBox::critical(this, tr("Error"), tr("Unable to laod any stacks. Please check parameters, and the stack file.") );
    project = old_project; return;
  }

  // NOTE: update_gl() populates vertex buffers and must run in a single (main) thread.
  if(conf.use_GL_buffers) project->update_gl(); 

  setWindowTitle("EDGE4D - " + QFileInfo(stackFile).fileName());

  view->setProject(project);
  view->setVol(project->stacks[0]);
  view->clearSelected();

  // Set slider positions.
  scrollTime->setRange(0, project->stacks.size() - 1); scrollTime->setSliderPosition(0);

  int width0 = project->stacks[0]->width(), height0 = project->stacks[0]->height(), depth0 = project->stacks[0]->depth();
  scrollXplane->setRange(0, width0 - 1);  scrollXplane->setSliderPosition(width0 / 2);
  scrollYplane->setRange(0, height0 - 1); scrollYplane->setSliderPosition(height0 / 2);
  scrollZplane->setRange(0, depth0 - 1);  scrollZplane->setSliderPosition(depth0 / 2);

  updateInfo();
  // Done setting slider positions.

  mTree->clear(); 
  mTree2Idx.clear();
  // TODO: Organize into noise, background, etc. at time point 0 get models 
  for(int idx = 0; idx < project->num_cells(0); idx++) {
    QString defaultText = tr("Model ") + QString::number(idx + 1);
    QTreeWidgetItem *modelItem = new QTreeWidgetItem(mTree, QStringList() << defaultText);
    modelItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    mTree2Idx[modelItem] = idx;
  }
  view->updateGL();

  // Display an information box with processing information.
  QString procTimeStr = "Run completed in " + QString::number(float(time.elapsed()) / 1e3) +  " seconds.";
  QMessageBox::information(this, tr("Processing time"), procTimeStr);
  if(old_project != NULL) delete old_project;
}

void EdgeWin::on_actionUpdate_Measurement_Parameters_triggered() {
  if(project == NULL) return;
  MeasureWin measure_win(project->conf, this); 
  if(measure_win.exec() == QDialog::Rejected) return; // User hit cancel.
}

MeasureWin::MeasureWin(edgeconfig &conf, QWidget *parent) : QDialog(parent), conf(conf) {
  setupUi(this);
  spinDFThreshold->setValue(conf.baz_thresh_pcent);
}

void MeasureWin::accept() {
  conf.baz_thresh_pcent = spinDFThreshold->value();
  // TODO Save .XML with measurement parameters.
  QDialog::accept();
}

// Application start point.
int main(int argc, char *argv[]) {
  // Launch main application.
  QApplication app(argc, argv);
  EdgeWin confocal; confocal.show();
  return app.exec();
}
