#ifndef __EDGE_HPP__
#define __EDGE_HPP__

#include "edgeconfig.hpp"
#include "mesh.hpp"
#include "volume.hpp"

#include <set>
#include <limits>
#include <fstream>

namespace edge {
  using namespace vol;

  struct dorsalfolds_t {
    dorsalfolds_t() {
      bazooka_pos = 0; 
      total_volume = 0;
      total_length = 0;
      apical_length = 0;
      total_surface_area = 0;
      above_baz_surface_area = 0;
      above_baz_volume = 0;
      cent_x = cent_y = cent_z = 0;
      basal_par_intensity = total_par_intensity = 0;
      bazooka_nonpatch_volume = bazooka_nonpatch_intensity = 0;
      bazooka_patch_volume = bazooka_patch_intensity = 0;
      basal_4voxel_par_volume = basal_4voxel_par_intensity = 0;
      basal_2voxel_par_volume = basal_2voxel_par_intensity = 0;
    }
    float cent_x, cent_y, cent_z;
    float apical_length, total_length;
    float total_surface_area, above_baz_surface_area; 
    float total_volume, above_baz_volume;
    float total_par_intensity, basal_par_intensity;
    float bazooka_pos;

    float bazooka_patch_volume, bazooka_patch_intensity;
    float bazooka_nonpatch_volume, bazooka_nonpatch_intensity;
    float basal_4voxel_par_volume, basal_4voxel_par_intensity;
    float basal_2voxel_par_volume, basal_2voxel_par_intensity;
  };

  // 3d meshes + meta data and measureemnts.
  struct image_cell {
    int idx;  // Index assigned to cell.
    bool displayed, selected;  // Flags controlling display/selected status.
    int trajectoryid; // Each cell mesh has an assigned trajectory ID (if < 0, no trajectory assigned)
    int bfs_label;

    edgeconfig &conf;  
    geom::mesh *mesh;

    vector<ivec3> voxels; // segmented voxels for a give cell
    vector<int> neighbors; // indicies in image_stack::cell_voxels of neighboring cells.
    vector<ivec3> surface_voxels; // surface voxel
    vector<int> surface_neighbors; // per surface voxel, indicies in image_stack::cell_voxels of neighboring cells
    
    // Organized by mesh face, indicies in image_stack::cell_voxels of neighboring cells.
    vector<int> mesh_face_neighbors;
    vector<float> pca_dims; // in microns PCA dimension

    // TODO: Figure out how this relates to VBOs because VBO objects
    // can have colors.
    vector<uint8> v_samples; // Sampled image volume on the mesh.
    int label_idx;  // index of label display 
    vector< vector<int> > labels; // Labels used for segmentation steps. 

    vec3 centroid, apical_pos, basal_pos;
    vector<vec3> avec;
    vector<vec3> v0, v1; // boundary edges, edges of trianges with no adjacencies

    plane3 p1, p2;

    // A mesh can have parts (individual meshes that are parts of the main mesh).
    vector<image_cell *> parts;  
    BVH2<image_cell> *parts_bvh;

    image_cell *nucleus;
    float vol_above_nuc, sa_above_nuc; // volume and surface area above nucleus.

    float surface_area, model_volume;     // Some measurements.
    float max_cell_bend;
    
    dorsalfolds_t dorsalfolds;
    float rel_angle, aspect_ratio, rel_aspect, anisotropy, anisotropy_AB;
    float depth_len, nuc_apical_dist;

    vector<float> neighbor_bend, neighbor_contact_sa;

    image_cell(int idx, geom::mesh *mesh, edgeconfig &c, vector<ivec3> &v, vector<int> &n, vector<ivec3> &sv, vector<int> &sn) : 
      conf(c), voxels(v), neighbors(n), surface_voxels(sv), surface_neighbors(sn)
    { 
      this->idx = idx; this->mesh = mesh; 
      bfs_label = -1;
      selected = displayed = false; 
      mesh->build_bvh();
      parts_bvh = NULL;  label_idx = 0;
      mesh->boundary_edges(v0,v1);
      surface_area = 0;
      model_volume = 0;
      rel_angle = 0;
      anisotropy = 0;
      anisotropy_AB = 0; 
      aspect_ratio = 1;
      depth_len = 0;
      trajectoryid = -1;
      centroid = mesh->centroid;
      nucleus = NULL;
      sa_above_nuc = vol_above_nuc = 0;
      max_cell_bend = 0;
      nuc_apical_dist = 0;
    }
    ~image_cell() { 
      if(mesh != NULL) delete mesh; 
      for(size_t p = 0; p < parts.size(); p++) delete parts[p];
      if(parts_bvh != NULL) delete parts_bvh;
    }

    // Change label index that is displayed. 
    void cycle_labels() { 
      if(labels.size() == 0) return; 
      label_idx = (label_idx + 1) % labels.size(); 
      mesh->boundary_edges(labels[label_idx], v0, v1);
    }

    void assign_neighbor_labels();

    void measure_dorsalfolds(vector<volume8 *> &vs, vector<volume8*> &dispvols);
    void measure_generic(vector<volume8*> &dispvols);

    void measure(vector<volume8 *> &vs, vector<volume8*> &dispvols) {
      switch(conf.analysis_id) {
      case analysis::DorsalFolds: measure_dorsalfolds(vs, dispvols); return;
      default: measure_generic(dispvols); return;
      }
    }
    void save_generic(string filename);

    void measure(image_cell *parent) { 
      surface_area = compute_sa();  // Compute surface area
      avec.clear();

      // Find index of face closest to mesh centroid.
      float mindsq = distance3sq(mesh->face_centroid(0), mesh->centroid);
      int fidx_center = 0;
      for(int fidx = 1; fidx < mesh->num_faces(); fidx++) {
	float dsq = distance3sq(mesh->face_centroid(fidx), mesh->centroid);
	if(dsq < mindsq) { mindsq = dsq; fidx_center = fidx; }
      }
      // Use "central face" to update mesh centroid.
      vec3 C = mesh->face_centroid(fidx_center); 
      mesh->centroid = C;
      mesh->pca3[0] = mesh->faces[fidx_center].normal;

      // Axial direction is the PCA of the tube (parent mesh).
      vec3 axial = parent->mesh->pca3[2];  

      vec3 proj_axial = mesh->pca3[2].proj(axial) + mesh->pca3[1].proj(axial);
      proj_axial.normalize();
      rel_angle = mesh->pca3[2].abs_angle(proj_axial) * (180.0 / M_PI);

      vector<vec3> edge; mesh->edge_vert(edge); if(edge.size() == 0) return;

      // Find the vector from the central face centroid to the shortest extent.
      vec3 e, emin; float mindot;
      emin = e = edge[0] - C;
      e.normalize();
      mindot = fabs(e.dot(mesh->pca3[2]));
      for(size_t i = 1; i < edge.size(); i++) {
	e = edge[i] - C;
	e.normalize();
	float dot = fabs(e.dot(mesh->pca3[2]));
	if(dot < mindot) { mindot = dot; emin = edge[i] - C; }
      }
      // avec.push_back(emin);

      // Use the shortest extent and the surface area to estimate a
      // long extent.  This assumes you can fit a ellipse
      // approximately to the data and assumes shortest extent is
      // accurately estimated.
      float short_extent = emin.length();
      float vox_surface_area = 0;
      vector<triangle> &tri = mesh->tribuf;    
      for(size_t fidx = 0; fidx < tri.size(); fidx++) {
	triangle t = tri[fidx]; vox_surface_area += t.area();
      }
      float long_extent_est = vox_surface_area / (M_PI * short_extent);
      aspect_ratio = long_extent_est / short_extent;

      vec3 emax; float maxdot;
      e = edge[0] - C;
      emin = e;
      emax = e;
      e.normalize();
      int idxmin = 0, idxmax = 0;
      maxdot = mindot = e.dot(proj_axial); 
      for(size_t i = 1; i < edge.size(); i++) {
	e = edge[i] - C;
	e.normalize();
	float dot = e.dot(proj_axial);
	if(dot < mindot) { mindot = dot; emin = edge[i] - C; idxmin = i; }
	if(dot > maxdot) { maxdot = dot; emax = edge[i] - C; idxmax = i; }
      }
      avec.push_back(emin); 
      avec.push_back(emax);
      float alpha = conf.voxel_alpha();
      vec3 amax = alpha * ( edge[idxmax] - C );
      vec3 amin = alpha * ( edge[idxmin] - C );
      anisotropy = amax.length() + amin.length();
    }

    void process_surface(volume8 &venh) {
      if(mesh == NULL) return;
      mesh->sample(v_samples, venh);  // Sample median filtered volume along tube surface.

     // Band pass that surface. 
      vector<int> bandpass1;
      const int surf_blurIter = 120;
      const int surf_threshold = 100;  
      mesh->surface_bandpass(bandpass1, v_samples, 0, 1, surf_threshold, surf_blurIter);
      labels.push_back(bandpass1);

      vector<int> bandpass2;  cout << "bandpass2" << endl;
      mesh->dilate_label(bandpass2, bandpass1, 3, 0, 1);
      labels.push_back(bandpass2);

      vector<int> seeds = bandpass2; cout << "seeds" << endl;
      int num_components = mesh->label_components(seeds);
      labels.push_back(seeds);
      
      vector<int> clean_seeds = seeds; cout << "clean_seeds" << endl;
      // NOTE: Clean seeds removes some labels. This will lead to gaps in the label set.     
      mesh->clean_small_comps(clean_seeds, num_components);
      labels.push_back(clean_seeds);

      vector<int> watershed = clean_seeds;  cout << "watershed" << endl;
      mesh->watershed(v_samples, watershed);    
      labels.push_back(watershed);

      // Set current label index and update boundary edges.
      label_idx = labels.size() - 1;
      mesh->boundary_edges(labels[label_idx], v0, v1);
    }

    // Use last label set to split cell shape into parts.
    void split_parts() {
      if(labels.size() == 0) return;
      vector<int> &labelvec = labels.back();
      // Get unique labels in label vector.
      set<int> ulabels;
      for(size_t i = 0; i < labelvec.size(); i++) ulabels.insert(labelvec[i]);
      // Collect face indices of each of these labels.
      vector<int> subidx;
      vector<vec3> vtable; vector<face> ftable;
      int part_id = 0;
      for(set<int>::iterator it = ulabels.begin(); it != ulabels.end(); it++) {
	// Get face indices of triangles with current label.
	subidx.clear();
	for(int fidx = 0; fidx < (int)labelvec.size(); fidx++) { if(labelvec[fidx] == *it) subidx.push_back(fidx); }
	// Build a mesh.
	mesh->as_mesh(subidx, vtable, ftable);
	vector<ivec3> vox_NA;
	vector<int> neighbor_NA;
	vector<ivec3> surface_voxels_NA; 
	vector<int> surface_neighbors_NA; 
	parts.push_back(new image_cell(part_id, new geom::mesh(vtable, ftable), conf, vox_NA, neighbor_NA, surface_voxels_NA, surface_neighbors_NA));
	part_id++;
	parts.back()->displayed = true;
      }
      cout << "parts.size()=" << parts.size() << endl;
      for(size_t i = 0; i < parts.size(); i++) {
	//parts[i]->measure(this);
	parts[i]->surface_area = parts[i]->compute_sa();
      }
      parts_bvh = new BVH2<image_cell>(parts);
    }

    void measure_lateral_bend();
    void measure_neighbor_contact_properties();

    image_cell *select_part(vec3 &o, vec3 &p) {
      if(parts.size() == 0 || parts_bvh == NULL) return NULL;
      // Find cells who's PARTS_BVH intersects the ray.
      vector<image_cell *> res; 
      parts_bvh->Query(res, o[0], o[1], o[2], p[0], p[1], p[2]); 
      if(res.size() == 0) return NULL;

      vector<image_cell *> isected; isected.reserve(res.size());
      vector<vec3>  isects; isects.reserve(res.size());
      for(size_t r = 0; r < res.size(); r++) {
	image_cell *cell = res[r];
	// If the cell has its displayed flag set, compute a ray
	// triangle intersection.
	if(cell->displayed) {
	  vec3 n;
	  if(cell->ray_intersect(n, o, p)) {
	    isected.push_back(cell);
	    isects.push_back(n);
	  }
	}
      }
      if(isects.size() == 0) return NULL;

      // Get the cell that has the closest intersection point to the ray.
      image_cell *icell = isected[0]; 
      vec3 isect = isects[0];
      float mindsq = distance3sq(isects[0], o);
      for(size_t i = 1; i < isected.size(); i++) {
	float dsq = distance3sq(isects[i], o);
	if(dsq < mindsq) {
	  mindsq = dsq;
	  icell = isected[i];
	  isect = isects[i];
	}
      }
      return icell;      
    }


    int analysis_id() { return conf.analysis_id; }
    int num_faces() { return mesh->num_faces(); }
    // n contains the intersection point if it exists.
    bool ray_intersect(vec3 &n, vec3 &o, vec3 &d) { int fidx; return mesh->ray_intersect(n, fidx, o, d); }

    // Cast a ray and get the face index of the intersected
    // triangle. Paint in a region around the intersected triangle.
    void label(vec3 &o, vec3 &d, int label) {
      if(labels.size() == 0) return;
      int fidx; vec3 n;
      // TODO: Get this connected with other label set.
      if(mesh->ray_intersect(n, fidx, o, d)) mesh->dfs_paint(labels.back(), fidx, label, 6);
    }
    float compute_sa() { 
      float surface_area = 0;
      float alpha = conf.voxel_alpha();
      vector<triangle> tri; mesh->as_triangles(tri);    
      for(size_t fidx = 0; fidx < tri.size(); fidx++) {
	triangle t = tri[fidx];
	// Map from isotropic voxel space to micron space.
	t[0] *= alpha; t[1] *= alpha; t[2] *= alpha;
	surface_area += t.area();
      }
      return surface_area;
    }
    vec3 get_centroid() { return conf.voxel_alpha() * centroid; }

    // Used by BVH2<> 
    float x() { return mesh->x(); } float y() { return mesh->y(); } float z() { return mesh->z(); }
    float min_x() { return mesh->min_x(); } float max_x() { return mesh->max_x(); }
    float min_y() { return mesh->min_y(); } float max_y() { return mesh->max_y(); }
    float min_z() { return mesh->min_z(); } float max_z() { return mesh->max_z(); }   
  };


  // Cell detection comparison.
  void static_analysis_internal(vector<int> &seg_status, 
				vector<int> &assigned_truth,
				int &missed_truth,
				ivec3 dims, 
				float truth_alpha, 
				float seg_alpha,
				vector< vector<ivec3> > &truth_voxels, vector< vector<ivec3> > &seg_voxels);

  enum neighbormode_t { Neighbor_All = 0, Neighbor_Xgreater = 1, Neighbor_Xless = 2 };

  struct image_stack {
    string basename;
    edgeconfig &conf;
    int frame_num;
    vector<volume16 *> vs16; // original 16-bit volume channels (freed during pre processing)
    vector<volume8 *> vs, dispvols;   // 8-bit isotropic 1:1:1 volumes for processing and 8-bit display volumes

    vector< ivec3 > edge_voxels;     // voxels from membrane processing
    vector< ivec3 > boundary_voxels; // voxels from outer boundary definition from repair_boundar()
    vector< ivec3 > mask_voxels;     // masked region voxels from repair_boundary()

    vector<image_cell *> cells, nuclei;  // cells/nucs segmented for display

    vector< vector<ivec3> > nuc_voxels; // fg = voxels of segmented nuclei (possibly merged)
    vector< vec3 > nuc_motion; // movement direction of nucleus
    vector<int> nuc_traj_ids; // trajectory IDs for nuclei
    map<int, int> nuc_traj_id2nuc; // map from trajectory id to nuc

    vector< vector<ivec3> > cell_voxels;  // Segmented cell voxels.

    vector< vector<int> > cell_neighbors; // Neighbor indices for each cell (index in cell_voxels).
    vector< vector<ivec3> > cell_surface_voxels; // Surface voxels of segmented cells.
    vector< vector<int> > cell_surface_neighbors; // Neighbors at each surface voxel ( index in cell_voxels ).

    vector< vector<ivec3> > cell_outline_voxels;

    vector< vec3 > cell_motion; // movement direction of EDT weighted cell centroid,  cellid -> motion vec
    vector<int> traj_ids; // Trajectory ID assignments to each cell voxel set.
    map<int, int> traj_id2cell; 
    vector< int > assigned_nucs; // Nucs assigned to each cell (index in nuc_voxels). 

    BVH2<image_cell> *bvh, *bvh_nuc; // BVH built on cells and nucs for clicking in GUI
    int iwidth, iheight, idepth, ichannels; // stack frame dimensions and channels

    int width() { return iwidth; } int height() { return iheight; }  int depth() { return idepth; }
    int nchannels() { return ichannels; }

    image_stack(string basename, vector<volume16 *> &vs16, edgeconfig &conf) : 
      basename(basename), conf(conf),  vs16(vs16) { 
      bvh_nuc = bvh = NULL; 
      ichannels = iwidth = iheight = idepth = 0;
      frame_num = -1;
    }
    ~image_stack() {
      // Free up memory used by volumes and cells. 
      for(size_t i = 0; i < vs16.size(); i++) delete vs16[i];
      for(size_t i = 0; i < vs.size(); i++) delete vs[i];
      for(size_t i = 0; i < dispvols.size(); i++) delete dispvols[i];
      for(size_t i = 0; i < cells.size(); i++) delete cells[i];      
      if(bvh != NULL) delete bvh; if(bvh_nuc != NULL) delete bvh_nuc;
    }

    void process() {
      pre_process();
      cout << "(image_stack) conf.analysis.id = "  << conf.analysis_id << endl;
      switch(conf.analysis_id) {
      case analysis::VolViewer: break;
	//case analysis::SPIMTest: process_nucs(); break;
      case analysis::CellShapes: process_edge(); break;
      case analysis::DorsalFolds: process_edge(); break;
      case analysis::NucsOnly: 
	if(ichannels != 1 && (conf.edge_channel - 1) >= 0) process_edge();
	process_nucs(); 
	break;	
      case analysis::NucsMembranes:  case analysis::StaticAnalysis:  
      case analysis::ManualGroundTruth:
	process_edge(); process_nucs(); break;
      default: break; 
      }           
      
    }

    volume8 *nuc_channel() {
      if(conf.keep_edge_channel_only) return NULL;  else return dispvols.at(conf.nuc_channel - 1);
    }
    volume8 *edge_channel() {
      if(conf.keep_edge_channel_only) return dispvols.at(0); else return dispvols.at(conf.edge_channel - 1);
    }

    // Process image stack, build 3d models of cells.
    void pre_process(); 
    volume8 *enhance(volume8 *v, int r1, float alpha1, int r2, float alpha2);
    void repair_boundary(volume8 &bvL, volume8 &bvH, boundaryconf_t &bconf);
    void build_meshes(vector< vector<ivec3> > &comps, vector< vector<int> > &neighbors, bool is_nuc = false);

    int clear_border() { return max(conf.r1,conf.r2) + 2*ceil(conf.sigma_low); }

    void compute_DOG_simple(volume8 *venh, volume8 *bv, int thresh, float sigma, bool keep_border = false);

    bool is_border_touch(vector<ivec3> &vox);

    void process_edge(); // Extracts membrane signal.
    void process_nucs(); // Extracts nuclei signal.
    void process_align(image_stack *stack);
    
    // Collect centroids of nuc voxels.
    void nuc_centroids(vector<vec3> &nuc_cent) {
      nuc_cent.resize(nuc_voxels.size());
      for(size_t n = 0; n < nuc_voxels.size(); n++) {
	nuc_cent[n] = geom::Centroid3(nuc_voxels[n]);
      }
    }


    // Computes centroids of cells, weighting based on EDT. Idea is that 
    // this favors better positioned centroids.
    void cell_centroids(vector<vec3> &cell_cent) {
      cell_cent.resize(cell_voxels.size());
      for(size_t c = 0; c < cell_voxels.size(); c++) {
	vector<ivec3> &vox = cell_voxels[c];
	// Convert voxels to vec3
	vector<vec3> vox_pts(vox.size());
	for(size_t i = 0; i < vox.size(); i++) vox_pts[i] = vox[i].to_vec3();
	cell_cent[c] = geom::Centroid3(vox_pts);

	vector<int> z_vals(vox.size());
	for(size_t i = 0; i < vox.size(); i++) z_vals[i] = vox[i].z();

	sort(z_vals.begin(), z_vals.end());
	int mid = z_vals.size() / 2;
	cell_cent[c].z() = z_vals[mid];
      }
    }

    image_cell *select_cell(vec3 &o, vec3 &p /*, image_cell *selected = NULL*/ ) {
      // TODO: Optionally handle selected cell (which is marked as undisplayed)
      // allow user to unselect a selected cell!!!!
      if(bvh == NULL) return NULL;

      // Find cells who's BVH intersects the ray.
      vector<image_cell *> res; 
      bvh->Query(res, o[0], o[1], o[2], p[0], p[1], p[2]); 
      if(res.size() == 0) return NULL;

      vector<image_cell *> isected; isected.reserve(res.size());
      vector<vec3>  isects; isects.reserve(res.size());
      for(size_t r = 0; r < res.size(); r++) {
	image_cell *cell = res[r];
	// If the cell has its displayed flag set, compute a ray
	// triangle intersection.
	if(cell->displayed) {
	  vec3 n;
	  if(cell->ray_intersect(n, o, p)) {
	    isected.push_back(cell);
	    isects.push_back(n);
	  }
	}
      }
      if(isects.size() == 0) return NULL;

      // Get the cell that has the closest intersection point to the ray.
      image_cell *icell = isected[0]; 
      vec3 isect = isects[0];
      float mindsq = distance3sq(isects[0], o);
      for(size_t i = 1; i < isected.size(); i++) {
	float dsq = distance3sq(isects[i], o);
	if(dsq < mindsq) {
	  mindsq = dsq;
	  icell = isected[i];
	  isect = isects[i];
	}
      }
      return icell;
    }

    void save_generic(string filename); 
    void save_analysis(string filename) {
      switch(conf.analysis_id) {
      case analysis::CellShapes: save_generic(filename); return;
      default: return;
      }
    }
    void save_raw_segmentation(string filename, bool clean_shapes = true);

    // Measure visble or selected cells.
    void measure_visible() {
      for(size_t c = 0; c < cells.size(); c++) { if(cells[c]->displayed) cells[c]->measure(vs, dispvols); }
    }
    void measure_selected() {
      for(size_t c = 0; c < cells.size(); c++) { if(cells[c]->selected) cells[c]->measure(vs, dispvols); }
    }

    void static_analysis_header(ofstream &output);
    void static_analysis(ofstream &output); 
    
    // Updates OpenGL vertex buffers. This MUST run in a (main) single thread.
    void update_gl() { 
      for(size_t c = 0; c < cells.size(); c++) cells[c]->mesh->update_gl(); 
      for(size_t n = 0; n < nuclei.size(); n++) nuclei[n]->mesh->update_gl(); 
    }

    void assign_nucs();
    void split_nucs(); 
    void find_cell_neighbors2();

    void ortho_comps_cells(vector< vector<ivec3> > &comp, int dim, int dim_pos);
    void ortho_comps_nucs(vector< vector<ivec3> > &comp, int dim, int dim_pos);

    void toggle_all_nucs() {
      for(size_t n = 0; n < nuclei.size(); n++) nuclei[n]->displayed = !nuclei[n]->displayed;
    }
    void match_nuclei_to_cells_displayed() {
      assign_nucs();
      for(size_t n = 0; n < nuclei.size(); n++) nuclei[n]->displayed = false;
      for(size_t a = 0; a < assigned_nucs.size(); a++) {
	if(assigned_nucs[a] >= 0 && cells[a]->displayed) {
	  nuclei.at(assigned_nucs[a])->displayed = true;
	}
      }
    }

    // Assigns nuclei trajectory IDs to cell shapes.
    void transfer_traj_ids() {
      assign_nucs();
      traj_ids.resize(assigned_nucs.size());
      for(size_t a = 0; a < assigned_nucs.size(); a++) {
	if(assigned_nucs[a] >= 0) traj_ids[a] = nuc_traj_ids.at(assigned_nucs[a]);
	else traj_ids[a] = -1;
      }
    }
    // Transfer nuc mesh data to cells. TODO: Clean up these interfaces.
    void transfer_nuc_meshes() {
      for(size_t a = 0; a < assigned_nucs.size(); a++) {
	if(assigned_nucs[a] >= 0) cells[a]->nucleus = nuclei.at(assigned_nucs[a]);
      }
    }

    void get_neighbors_traj(vector<int> &neighbor_traj, neighbormode_t neighbormode, vector<int> &order, vector<int> &cur_traj);
    void get_neighbors_traj(vector<int> &neighbor_traj, int order, vector<int> &cur_traj) {
      vector<int> order_vec; order_vec.push_back(order);
      get_neighbors_traj(neighbor_traj, Neighbor_All, order_vec, cur_traj);
    }
    void get_neighbors(vector<int> &neighbor_idx, neighbormode_t neighbor_mode, vector<int> &order, vector<int> &cur_idx);
  };

  struct hyperstack {
    string basename;
    edgeconfig conf;  // Global configuration information for hyperstack.
    vector<image_stack *> stacks; // Multiple stacks in a hyperstack.
    vector<string> files;

    vector<int> nuc_traj_lengths, traj_lengths; // id -> trajectory length  mapping.

    hyperstack(string filename, edgeconfig &conf) : conf(conf) {
      files.push_back(filename);
    }
    hyperstack(vector<string> &filelist, edgeconfig &conf) : conf(conf) {
      files = filelist;
    }
    
    // Free up memory associated with stacks.
    ~hyperstack() { for(size_t i = 0; i < stacks.size(); i++) delete stacks[i];  }
    void run_segmentation(bool split_nucs = false);
    void track_nucs_centroid2();
    void track_cells_centroid2();
    void load_tiff_stacks();
    void process();

    int num_cells(int stack) { return stacks.at(stack)->cells.size();  }
    int num_stacks() { return stacks.size(); }
    // Measures selected or visible meshes given stack index.
    void measure_visible(int stack_idx) { stacks.at(stack_idx)->measure_visible(); }
    void measure_selected(int stack_idx) { stacks.at(stack_idx)->measure_selected(); }

    void update_gl() { for(size_t t = 0; t < stacks.size(); t++) stacks[t]->update_gl();   }

    void static_analysis(string filename) {
      if(stacks.size() == 0) return;
      ofstream output(filename.c_str());
      output << "frame\t" << "time\t";
      stacks[0]->static_analysis_header(output);
      output << endl;
      // Output a line for each stack.
      for(size_t t = 0; t < stacks.size(); t++) {
	output << t + conf.startframe << '\t';
	output << (t + conf.startframe) * conf.time_step << '\t';
	stacks[t]->static_analysis(output);
	output << endl;
      }
    }

    void global_analysis(string filename) {
      if(conf.analysis_id == analysis::StaticAnalysis) static_analysis(filename);
    }

    void save_analyze_trajectories(string filename);
    void save_neighbor_analysis(string filename);
    void save_contact_analysis(string filename);

    void save_neighbor_swaps(string filename);

    void save_trajectory_lengths(string filename) {
      ofstream output(filename.c_str());
      output << "id" << '\t' << "traj.length" << endl;
      for(size_t i = 0; i < traj_lengths.size(); i++) {
	if(traj_lengths[i] >= 2) {  output << i << '\t' <<  traj_lengths[i] << endl; }
      }
    }

    void save_stable_neighbors(string filename);
    void save_neighbor_counts(string filename);

    void measure_trajectories_threaded();
    void save_dorsalfolds(string filename, vector<int> &traj_ids);
    void save_dorsalfolds(string filename) {
      vector<int> traj_ids_all;
      for(int traj_id = 0; traj_id < (int)traj_lengths.size(); traj_id++) {
	if(traj_lengths[traj_id] >= 1) traj_ids_all.push_back(traj_id);
      }
      save_dorsalfolds(filename, traj_ids_all);
    }

    void save_nuc_cents(const string &filename);
  };
};


#endif // __EDGE_HPP__
