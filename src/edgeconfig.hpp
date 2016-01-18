#ifndef __EDGECONFIG_HPP__
#define __EDGECONFIG_HPP__

#include <vector>
#include <string>

using namespace std;

// TODO: Need to somehow separate GUI parameters and state information
// from data processing parameters.

namespace analysis {
  // Processing pipelines.
  const int VolViewer = 0; 
  const int CellShapes = 1; 
  const int DorsalFolds = 2;
  const int NucsOnly = 3;
  const int NucsMembranes = 4;

  const int StaticAnalysis = 5; // get through these analyses one by one
  const int ManualGroundTruth = 6;
  const int ACME_Compare = 7;

  // Static analysis error types. 
  const int StaticNoError = 0;
  const int StaticJunk = 1;
  const int StaticSplit = 2;
  const int StaticMerge = 3;
  const int StaticSplitMerge = 4;
};

struct boundaryconf_t {
  int startframe, endframe;
  int boundary_mincomp;
  float boundary_fmin;
  bool boundary_all;
  int thresh_high;
  float sigma_high;
  int r_repair, r2_repair; // radius used for repairing outer boundaries
  boundaryconf_t() {
    startframe = endframe = 0;
    boundary_mincomp = 200;
    boundary_fmin = 0.7;
    boundary_all = false;
    thresh_high = 5; 
    sigma_high = 2.5;
    r2_repair = r_repair = 15;
  }
};

struct edgeconfig {
  bool use_GL_buffers;

  // System settings
  int threads, analysis_id, hs_threads; // CPU threads and type of analysis conducted on data
  bool run_refinement, run_mesh;
  bool keep_processing;
  bool keep_edge_channel_only;

  int edge_channel; // channel 1, 2, or 3 contains the edge (1-indexed!! so subtract 1)
  int nuc_channel; // channel is 1-indexed

  // Hyperstack metadata loaded from tiff file.
  int frames, slices, channels;
  float time_step; // (usually in seconds?)
  bool override_tiff;   // If true, user TIFF parameters should override hyperstack metadata.
  // User configured parameters for hyperstacks (movies).
  int startframe, endframe, stepframe;

  // Geometry and volume scaling information
  float voxelXY, voxelX, voxelY, voxelZ, scalef; 
  bool scaleToZ; // Adjusts x,y dimension to match z dimension.

  // Membrane enhancment parameters.
  int hist_nx, hist_ny, hist_nz;  // Number of histograms in volume.
  float fmin, fmax; // Histogram adjustment parameter.
  int r1, r2;  // rank filter radius
  float alpha1, alpha2; // alpha = fraction used in rank filter enhancement 
  int r1_nuc, r2_nuc;
  float alpha1_nuc, alpha2_nuc;
  int r_nonedge; // radius of median filter applied to non-edge volumes

  vector<boundaryconf_t> boundaryconf;

  // DOG/bandpass computation parameters.
  float sigma_low, sigma_nuc;
  int thresh_low, thresh_nuc;

  int nuc_split_dist;
  // min/max volume of nucleus in voxels, TODO: should display these parameters.
  int min_nuc_vol, max_nuc_vol; 
  float nearestk_dist;
  bool repair_nucs;

  // Segmentation parameters.
  int max_hole_rad, min_comp_vol, max_comp_vol, internal_limit, noise_comp_vol;
  float neighbor_alpha;

  // Mesh construction and refinement parameters.
  int mcdim;
  float refine_sigma, refine_dist, refine_alpha, refine_stepsize;
  int refine_iter;
  // -----------------------------

  // NM2010 params
  bool run_NM2010;
  bool use_ASF; int maxASF_rad;
  float filter_sigma;
  int bgseed_thresh;
  int h_seed;
  int error_vox_thresh;

  // Measurement parmeters
  float baz_thresh_pcent;

  float ACME_average_radius;
  float ACME_plate_measure;
  float ACME_voting_field;
  float ACME_membrane_segmentation;

  // Constructor sets default paramters.
  edgeconfig() {
    // System parameter defaults.
    use_GL_buffers = true;
    hs_threads = 4;
    threads = 4;
    analysis_id = 0;
    edge_channel = 1;
    nuc_channel = 2;
    run_refinement = true;
    run_mesh = true; 
    keep_edge_channel_only = false;
    keep_processing = false;

    // TIFF parameter defaults.
    frames = 1; slices = 0; channels = 1;     
    override_tiff = false;
    startframe = endframe = 0; 
    stepframe = 1; 
    voxelX = voxelY = voxelXY = voxelZ = scalef = 1.0;
    scaleToZ = true;
    time_step = 1;

    // Membrane enhancement parameters.
    hist_nx = hist_ny = 1; hist_nz = 8; 
    fmin = fmax = 0;
    r1 = 4; r2 = 4;
    alpha1 = 0.95; alpha2 = 0.05;

    // Nucleus enhancment parameters.
    r1_nuc = 3; r2_nuc = 3;
    alpha1_nuc = 0.05; alpha2_nuc = 0.95;
    repair_nucs = true;

    r_nonedge = 0;

    boundaryconf_t global_boundaryconf;
    boundaryconf.push_back(global_boundaryconf);


    sigma_low = 1.25;  sigma_nuc = 2.25;
    thresh_low = 2;    thresh_nuc = 8;

    min_nuc_vol = 200; max_nuc_vol = 2000;

    neighbor_alpha = 0.1;

    nuc_split_dist = 7; // Search distance for hausdorff matching, show parameter in GUI
    nearestk_dist = 16;

    max_hole_rad = 8; 
    min_comp_vol = 1000; 
    max_comp_vol = 20000; 
    noise_comp_vol = 100;
    internal_limit = 600;

    mcdim = 2; 
    refine_sigma = 3; refine_dist = 2; 
    refine_iter = 10;
    refine_alpha = 0.8; refine_stepsize = 0.2;

    run_NM2010 = false; 
    use_ASF = false; 
    maxASF_rad = 2;
    filter_sigma = 1.0;
    bgseed_thresh = 4;
    h_seed = 4; // h_seed should be less than equal to bgseed_thresh?
    error_vox_thresh = 100;

    ACME_average_radius = 0.3; 
    ACME_plate_measure = 0.7;
    ACME_voting_field = 1.0;
    ACME_membrane_segmentation = 1.0;

    // measurement parameters
    baz_thresh_pcent = 0.99;
  }

  // Computes the scale factor using the above information.  This
  // scale factor assures that a 1:1:1 volume ratio.
  void get_scale(float &sx, float &sy, float &sz){
    if(scaleToZ) { sx = sy = scalef * (voxelXY / voxelZ); sz = scalef; }
    else { sx = sy = scalef; sz = scalef * (voxelZ / voxelXY);  }
  }

  // Get multiplicative factor for converting voxel units to microns. 
  float voxel_alpha() { 
    float alpha = 1; 
    if(scaleToZ) alpha = voxelZ; else alpha = voxelXY;
    alpha /= scalef;
    return alpha;
  }

  // Loads metadata from a TIFF stack file created with ImageJ. 
  bool load_tiff(string tiffstack);

  bool load_xml(string filename);
  bool save_xml(string filename);

  void getBoundaryConfig(boundaryconf_t &bconf, int frame_num) {
    bconf = boundaryconf.at(0);
    for(size_t b = 0; b < boundaryconf.size(); b++) {
      if(boundaryconf[b].startframe <= frame_num && frame_num <= boundaryconf[b].endframe) {
	bconf = boundaryconf[b];
      }
    }
  }

};

#endif // __EDGECONFIG_HPP__
