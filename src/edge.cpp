#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <set>
#include <queue>
#include <sstream>

#include "marchingcubes.hpp"
#include "edge.hpp"
#include "procthreads.hpp"

#include <QThread>
#include <QFileInfo>

// TODO: get rid of vec3 
#include <QMatrix4x4> 

#include "util.hpp"

using namespace vol;
using namespace util;


namespace edge {
  void image_stack::pre_process() {
    if(vs16.size() == 0) { cout << "no channels to pre-process" << endl; return; }
    vs.resize(vs16.size());

    // TODO: For 5d data, we should provide a consistent stretch to each volume.
    for(size_t channel = 0; channel < vs.size(); channel++) {
      // Stretch 16-bit values in entire volume so they span whole range [0,65535] range.
      // Don't do this if the original was 8-bit (at load I assume everything is 16-bit). 
      // NOTE: We may want to use a consistent operation here for all the channels and/or volumes
      // especially when intensity measurements are needed.
      if(vs16[channel]->is8bitOriginal == false) stretch16(vs16[channel], conf.threads);

      // Now convert to 8-bit.
      vs[channel] = new volume8(*vs16[channel]);
      delete vs16[channel]; // Free up memory used by 16-bit volumes.
    }
    vs16.clear();  // Clear vector.

    // Pre-proceses all channels. 
    for(size_t c = 0; c < vs.size(); c++) {
      // Scale the volume so voxels are all the same size. Scale to
      // match the z-axis by stretching/compressing the x and y axes
      // or to match the x and y axes by stretching/compressing the
      // z-axis.  
      if(conf.voxelX != conf.voxelY) {
	if(conf.voxelX < conf.voxelY) vs[c] = scaleUpXY(vs[c], 1.0, conf.voxelY / conf.voxelX, conf.threads);
	else vs[c] = scaleUpXY(vs[c], conf.voxelX / conf.voxelY, 1.0, conf.threads);
      }

      // Scale volume so isotropic.
      if(conf.scaleToZ) vs[c] = scaleToZ(vs[c], conf.voxelXY / conf.voxelZ, conf.threads); 
      else vs[c] = scaleToXY(vs[c], conf.voxelZ / conf.voxelXY, conf.threads);

      // Scale volume globally. 
      if(conf.scalef != 1.0f) vs[c] = scaleXYZ(vs[c], conf.scalef, conf.threads);

      /*if(frame_num == 1 || frame_num == 3) {
	  volume8 *src = vs[c];
	  volume8 *dst = new volume8(src->width, src->height, src->depth);
	  dst->mirrorx(*src); // mirror x-axis to align two views
	  delete src;
	  vs[c] = dst;
      }
      if(frame_num == 2 || frame_num == 3) {
	cout << "rotate90" << endl;
	volume8 *src = vs[c];
	volume8 *dst = new volume8(src->width, src->height, src->depth);
	dst->rotate90minus(*src);  // -90 rotatation to align 
	delete src;
	vs[c] = dst;
	}*/
      
      cout << "channel = " << c << endl;
      cout << "scaled dimensions:" << vs[c]->width << " x " << vs[c]->height << " x " << vs[c]->depth << endl;
      cout << "voxel_alpha() = " << conf.voxel_alpha() << endl;
      float alpha = conf.voxel_alpha();
      cout << "scale sized in um: " << float(vs[c]->width) * alpha << " x " << float(vs[c]->height) * alpha << " x " << float(vs[c]->depth) * alpha << endl;

      if(!conf.keep_edge_channel_only) {
	dispvols.push_back(new volume8(*vs[c]));    // Keep 8-bit scaled volume for display.
	// Enhance display only when fmax is applied.
	//if((conf.fmin > 0 && conf.fmax > 0) || (conf.fmin == 0 && conf.fmax > 0)) 
	//happly(dispvols.back(), conf.hist_nx, conf.hist_ny, conf.hist_nz, conf.fmin, conf.fmax, conf.threads);
      }

    }
    // Keep edge channel only, if requested.
    if(conf.keep_edge_channel_only) dispvols.push_back(new volume8(*vs.at(conf.edge_channel-1)));

    // Save dimensions of image stack using width height and depth of 0th channel.
    // TODO: Need to check all volumes are the same.
    iwidth = vs[0]->width; iheight = vs[0]->height; idepth = vs[0]->depth;
    ichannels = vs.size();

    // Keep originals in a single channel if requested.
    //if(ichannels == 2 && conf.keep_processing) {
    //dispvols.push_back(new volume8(*vs[0])); dispvols.push_back(new volume8(*vs[1]));
    //}

    // Apply median filter to non-edge channel, if requested.
    if(conf.r_nonedge > 0) {
      cout << "r_nonedge=" << conf.r_nonedge << endl;
      vector<volume8 *> non_edge_new;
      for(size_t c = 0; c < vs.size(); c++) {
	// Rank filter non-edge volumes if requested (r_nonedge > 0).

	// TODO: Maybe this is better here since I can measure one cell at a time.
	if((int)c != (conf.edge_channel-1)) {
	  non_edge_new.push_back(new volume8(*vs[c]));
	  rankfilt(*(non_edge_new.back()), conf.r_nonedge, 0.5, conf.threads);
	}
      }
      for(size_t n = 0; n < non_edge_new.size(); n++) {
	vs.push_back(non_edge_new[n]);
	if(conf.keep_edge_channel_only == false && conf.keep_processing && conf.r_nonedge > 0)  
	  dispvols.push_back(new volume8(*non_edge_new.at(n)));
      }
    }
  }

  // Note: Caller must free memory allocated by returned volume.
  // Note: It is better to use a small r1 and r2 as possible. alpha1 = 0.95 and alpha1 = 0.05
  //       seem to work best for membranes.
  volume8 *image_stack::enhance(volume8 *v, int r1, float alpha1, int r2, float alpha2) {
    volume8 *venh = new volume8(*v); 

    // Apply local histogram adjustment. 
    if(conf.fmin > 0 && conf.fmax == 0) 
      happly(venh, conf.hist_nx, conf.hist_ny, conf.hist_nz, conf.fmin, 1e-6, conf.threads);
    else if(conf.fmin > 0 && conf.fmax > 0) 
      happly(venh, conf.hist_nx, conf.hist_ny, conf.hist_nz, conf.fmin, conf.fmax, conf.threads);

    // Perform rank filtering to remove high frequency detail while also preserving edges.
    cout << "rank filtering" << endl;
    if(conf.keep_processing && conf.fmin > 0) { dispvols.push_back(new volume8(*venh)); }
    if(conf.r1 > 0 && conf.r2 > 0) {
      rankfilt(*venh, r1, alpha1, conf.threads); // Rank filter to enhance edges.
      rankfilt(*venh, r2, alpha2, conf.threads); // Rank filter again to thin enhanced edges.
    }
    // Keep displayable version of enhanced volume, if requested.
    if(conf.keep_processing) dispvols.push_back(new volume8(*venh));     

    return venh;
  }

  // Simple DOG using low threshold and sigma_low.
  void image_stack::compute_DOG_simple(volume8 *venh, volume8 *bv, int thresh, float sigma, bool keep_border) {
    volume8 *v2 = new volume8(*venh);

    for(int dim = 0; dim < 3; dim++) blur(v2, (gauss_dim)dim, sigma, conf.threads);
    bv->subtract_and_clip(*venh, *v2, thresh);    
    delete v2;

    if(keep_border == false) bv->clear_border(clear_border());
  }


  // Use the high threshold DOG to find missing boundaries. Write
  // these boundaries to a low threshold DOG.

  // Only applies boundaries above and below the tissue structure.
  void image_stack::repair_boundary(volume8 &bvL, volume8 &bvH, boundaryconf_t &bconf) {
    int r = max(1, bconf.r_repair), r2 = max(1, bconf.r2_repair);  // Get repair radius.
    int rsq = r*r, r2sq = r2*r2; // Get repair radius squared.

    if(bconf.boundary_all) {
      bvL.pad(r+1, 0); bvH.pad(r+1, 0); // Pad volumes
    }
    else {
      bvL.pad_AB(r+1, 0); bvH.pad_AB(r+1, 0); // Pad volumes above and below.
    }

    // Compute EDT using high threshold DOG volume.
    volume32 *EDT = new volume32(bvH.width, bvH.height, bvH.depth);
    compute_edt(*EDT, bvH, 255, conf.threads);

    // Threshold EDT.  Look at "boundary" voxels.  Find regions that
    // are far from edge information.
    volume8 front(bvH.width, bvH.height, bvH.depth);
    EDT->above_threshold(front, rsq);

    // Compute the EDT from those far regions.
    compute_edt(*EDT, front, 255, conf.threads);

    // Get region rsq back from those far regions. 
    EDT->below_threshold(front, r2sq);
    delete EDT;

    for(int z = 0; z < front.depth; z++) 
      for(int y = 0; y < front.height; y++)
	for(int x = 0; x < front.width; x++) { if(front(x,y,z) == 255) bvL(x,y,z) = 0; }

    // Add a boundary.
    boundary_voxels.clear();
    vector<ivec3> adj(6);
    for(int z = 1; z < front.depth - 1; z++) 
      for(int y = 1; y < front.height - 1; y++)
	for(int x = 1; x < front.width - 1; x++) {
	  if(front(x,y,z) == 0) {
	    front.nhood6_unsafe(adj, ivec3(x,y,z));
	    int fgsum = 0;
	    for(size_t i = 0; i < adj.size(); i++) if(front(adj[i]) == 255) fgsum++;
	    if(fgsum > 0) {
	      bvL(x,y,z) = 255;
	      ivec3 b(x,y,z);
	      if(bconf.boundary_all) { ivec3 br(r+1,r+1,r+1); boundary_voxels.push_back(b - br); }
	      else {  // above/below (AB) padding
		ivec3 br(0,0,r+1);boundary_voxels.push_back(b - br);
	      }
	    }
	  }
	}


    // Remove padding.
    if(bconf.boundary_all) { bvL.unpad(r+1); bvH.unpad(r+1); front.unpad(r+1); }
    else {  // Remove above and below padding.
      bvL.unpad_AB(r+1); bvH.unpad_AB(r+1); front.unpad_AB(r+1);  
    }

    // Collect mask voxels in boundary region as well as
    if(bconf.boundary_all) {
    for(int z = 0; z < front.depth; z++) 
      for(int y = 0; y < front.height; y++)
	for(int x = 0; x < front.width; x++) { if(front(x,y,z) == 255) mask_voxels.push_back(ivec3(x,y,z)); }
    }
    else {
      // Add border region to mask voxels.
      int clear_brdr = clear_border();
      for(int z = 0; z < front.depth; z++) 
	for(int y = 0; y < front.height; y++)
	  for(int x = 0; x < front.width; x++) { 
	    if(front(x,y,z) == 255 || 
	       x < clear_brdr || y < clear_brdr  || z < clear_brdr ||
	       x >= front.width - clear_brdr || y >= front.height - clear_brdr || z >= front.depth  - clear_brdr ) {
	      mask_voxels.push_back(ivec3(x,y,z)); 
	      bvL(x,y,z) = 0;
	    }
	  }
      // Apply the border region.
      for(int z = 0; z < bvL.depth; z++) 
	for(int y = 0; y < bvL.height; y++)
	  for(int x = 0; x < bvL.width; x++) { 
	    if(front(x,y,z) != 255) {
	      if(x == clear_brdr || y == clear_brdr  || z == clear_brdr ||
		 x == front.width - clear_brdr - 1 || y == front.height - clear_brdr - 1  || z == front.depth  - clear_brdr - 1 ) {
		bvL(x,y,z) = 255;
	      }
	    }
	  }
    }

    // Save boundary repaired volume.
    if(conf.keep_processing) dispvols.push_back(new volume8(bvL));
  }


  
  // Use EDT + fractional Hausdorff distance. EDT values sorted with distance at f = 0.
  // Try several positions of previous voxelized shape.
  void image_stack::process_nucs() {
    // Enhance nuclei volume. 
    int w = iwidth, h = iheight, d = idepth;
    volume8 *venh = enhance(vs.at(conf.nuc_channel - 1), conf.r1_nuc, conf.alpha1_nuc, conf.r2_nuc, conf.alpha2_nuc); 

    volume8 *bv = new volume8(w,h,d);
    // Compute a simple DOG using a largeish sigma, and high threshold to find nuc regions.
    compute_DOG_simple(venh, bv, conf.thresh_nuc, conf.sigma_nuc, true);
    if(conf.keep_processing) dispvols.push_back(new volume8(*venh)); 
    delete venh;

    // Clear edge voxels to break some merged nucs. Also keep
    // edge voxels to allow reconstruction during temporal coherence processing.
    for(size_t e = 0; e < edge_voxels.size(); e++) bv->v(edge_voxels[e]) = 0;
    // Also clear masked region in nuc image.
    for(size_t m = 0; m < mask_voxels.size(); m++) bv->v(mask_voxels[m]) = 0;
    
    if(conf.keep_processing) dispvols.push_back(new volume8(*bv));

    // Threshold volume of nuc regions.
    bv->components(nuc_voxels, conf.min_nuc_vol, 255); 
    
    // Apply max volume threshold if analyzing nucs only.
    if((conf.analysis_id == analysis::NucsOnly ||
	conf.analysis_id == analysis::SPIMTest)
       && conf.repair_nucs == false) {
      vector< vector<ivec3> > nuc_voxels2;
      for(size_t n = 0; n < nuc_voxels.size(); n++) { 
	if((int)nuc_voxels[n].size() >= conf.min_nuc_vol && 
	   (int)nuc_voxels[n].size() <= conf.max_nuc_vol) nuc_voxels2.push_back(nuc_voxels[n]);
      }
      nuc_voxels = nuc_voxels2;
    }

    // Save display of nuc and membrane data.
    if(conf.keep_processing) {
      bv->fill(0);
      for(size_t e = 0; e < edge_voxels.size(); e++) bv->v(edge_voxels[e]) = 255;
      for(size_t n = 0; n < nuc_voxels.size(); n++) {
	vector<ivec3> &nuc = nuc_voxels[n];
	for(size_t i = 0; i < nuc.size(); i++) bv->v(nuc[i]) = 165;
	// Also draw outer boundary.
	for(size_t i = 0; i < boundary_voxels.size(); i++) bv->v(boundary_voxels[i]) = 130;
      }
      dispvols.push_back(new volume8(*bv));
    }

    delete bv; 
  }

  void image_stack::process_align(image_stack *stack2) {
    if(vs.size() == 0 || stack2->vs.size() == 0) return; // no channels?

    volume8 transformed(stack2->vs[0]->width, stack2->vs[0]->height, stack2->vs[0]->depth);

    QMatrix4x4 bestT;
    double bestcorr = 0;
    //    for(float dx = -6; dx <= 6; dx += 1) {
    //      for(float dy = -6; dy <= 6; dy += 1) {
    //    	for(float dz = -6; dz <= 6; dz += 1) {
    float dx = -6, dy = -1,  dz = 2;
    	  cout << dx << "\t" << dy << "\t" << dz << endl;
	  for(float a = -8; a <= 8; a += 2) {
	    for(float b = -8; b <= 8; b += 2) {
	      for(float c = -8; c <= 8; c += 2) {
		QMatrix4x4 T;
		T.translate(dx,dy,dz);
		QQuaternion rot = QQuaternion::fromEulerAngles(a,b,c);
		T.rotate(rot);
		transformed.fill(0);
		transformed.ridged(*stack2->vs[0], T);
		double corr = vs[0]->correlation(transformed);
		if(corr > bestcorr) {
		  cout << corr << endl;
		  bestcorr = corr;
		  bestT = T;
		}
	      }
	    }
	    //	  }
	  //	}
  //      }
    }
    
    qWarning() << bestT << endl;
      // Keep transformed volume.
    transformed.fill(0);
    transformed.ridged(*stack2->vs[0], bestT);
      
    dispvols.push_back(new volume8(transformed));
  }
  
  void image_stack::process_edge() {
    int w = iwidth, h = iheight, d = idepth;
    // Enhance by rank filtering.
    volume8 *venh = enhance(vs.at(conf.edge_channel - 1), conf.r1, conf.alpha1, conf.r2, conf.alpha2); 
    volume8 *bvL = new volume8(w,h,d); 
    compute_DOG_simple(venh, bvL, conf.thresh_low, conf.sigma_low);   

    
    volume8 *bvH = new volume8(w,h,d); 
    venh->copy(*vs.at(conf.edge_channel - 1));

    boundaryconf_t bconf;
    conf.getBoundaryConfig(bconf, frame_num);

    // Apply histogram cutoff to found boundary enclosing regions.
    happly(venh, 1, 1, 1, bconf.boundary_fmin, 1e-6, conf.threads);
    if(conf.keep_processing) dispvols.push_back(new volume8(*venh));

    compute_DOG_simple(venh, bvH, bconf.thresh_high, bconf.sigma_high);  
    delete venh; venh = NULL;

    // Thin low threhsold DOG.
    if(conf.keep_processing) dispvols.push_back(new volume8(*bvL));
    thin(*bvL, thinSURFACE, conf.threads);
    if(conf.keep_processing) dispvols.push_back(new volume8(*bvL));

    if(conf.keep_processing) dispvols.push_back(new volume8(*bvH));

    // Clear small components in boundary fg region.
    vector< vector<ivec3> > compsH;
    bvH->components(compsH, 0, 255, true);
    for(size_t i = 0; i < compsH.size(); i++){
      vector<ivec3> &comp = compsH[i];
      if((int)comp.size() <= bconf.boundary_mincomp) {  
	for(size_t c = 0; c < comp.size(); c++)  bvH->v(comp[c]) = 0; 
      }
    }

    if(conf.keep_processing) { // Keep DOG data for display if requested.
      dispvols.push_back(new volume8(*bvH)); dispvols.push_back(new volume8(*bvL)); 
    }

    // This step and those above, try to get rid of background noise.
    repair_boundary(*bvL, *bvH, bconf);
    delete bvH; 

    // Remove small thinned components ----------------------------
    vector< vector<ivec3> > thincomps;
    bvL->components(thincomps, 0, 255, true);
    for(size_t i = 0; i < thincomps.size(); i++){
      vector<ivec3> &comp = thincomps[i];
      if((int)comp.size() <= conf.noise_comp_vol) { for(size_t c = 0; c < comp.size(); c++)  bvL->v(comp[c]) = 0; }
    }
    if(conf.keep_processing) dispvols.push_back(new volume8(*bvL));
    // -------------------------------------------------------------

    // Save fg, edge voxels.
    for(int z = 0; z < d; z++)
      for(int y = 0; y < h; y++) 
	for(int x = 0; x < w; x++) { if(bvL->v(x,y,z) == 255) edge_voxels.push_back(ivec3(x,y,z)); }

    delete bvL;
  }

  // Given a set of components build and refine meshes. 
  void image_stack::build_meshes(vector< vector<ivec3> > &comps, vector< vector<int> > &neighbors, bool is_nuc) {
    cout << "build_meshes(): comps.size()=" << comps.size() << endl;

    // 10. Triangulate remaining components. TODO: Make this multithreaded.
    vector<mesh *> meshes;
    for(size_t c = 0; c < comps.size(); c++) {
      // Create a small binary volume for each component in the larger binary volume.
      volume8 bvtest(comps[c], conf.mcdim + 1 + 4);
      vector< vector<ivec3> > holes;
      // Get rid of any small background regions.
      bvtest.components(holes, 0, 0);
      sort(holes.begin(), holes.end(), cmp_size_gt<ivec3>());
      for(size_t h = 1; h < holes.size(); h++) { // skip largest, since it is background region
	for(size_t i = 0; i < holes[h].size(); i++) bvtest(holes[h][i]) = 255;
      }

      // Triangulate the volume.
      vector<vec3> vtable; vector<face> ftable;

      // NOTE: The 127.5 is critical to put the verticies right in the middle of the edge.
      marching::MarchingCubes(127.5, bvtest, vtable, ftable, conf.mcdim);

      vec3 shift(bvtest.x0, bvtest.y0, bvtest.z0);
      for(size_t vidx = 0; vidx < vtable.size(); vidx++) vtable[vidx] += shift;

      meshes.push_back(new mesh(vtable, ftable)); 
    }

    // 11. Refine triangulation to match image intensity.
    if(conf.run_refinement) {
      // Shit this used  to be 0!!!! Need to make sure edge channel is used.
      volume8 vblur = *vs.at(conf.edge_channel-1); 
      for(int dim = 0; dim < 3; dim++) blur(&vblur, (gauss_dim)dim, px2sigma(conf.refine_sigma), conf.threads);
      if(conf.keep_processing) dispvols.push_back(new volume8(vblur));
      vector<mesh *> refine_meshes;
      for(size_t i = 0; i < meshes.size(); i++) refine_meshes.push_back(meshes[i]);
      if(conf.analysis_id == analysis::NucsOnly || is_nuc){
	for(size_t r = 0; r < refine_meshes.size(); r++) refine_meshes[r]->smooth(conf.refine_iter);
      }
      else {
	refine(&vblur, 
	       conf.refine_dist, conf.refine_iter, conf.refine_alpha, conf.refine_stepsize, 
	       refine_meshes, conf.threads);
      }
    }

    // Set mesh data of individual cells.
    if(is_nuc) {
      for(size_t idx = 0; idx < meshes.size(); idx++) {
	// TODO: These NA entries are a hack. Clean up at some point.
	vector<int> neighbors_NA;
	vector<ivec3> surface_voxels_NA; 
	vector<int> surface_neighbors_NA; 
	nuclei.push_back(new image_cell(idx, meshes[idx], conf, comps.at(idx), neighbors_NA, surface_voxels_NA, surface_neighbors_NA));
	if(nuc_traj_ids.size() > 0) nuclei.back()->trajectoryid = nuc_traj_ids.at(idx);
      }      
    }
    else {
      for(size_t idx = 0; idx < meshes.size(); idx++) {
	if(neighbors.size() == 0) {
	  // 
	  vector<int> neighbors_NA;
	  vector<ivec3> surface_voxels_NA; 
	  vector<int> surface_neighbors_NA; 
	  cells.push_back(new image_cell(idx, meshes[idx], conf, comps.at(idx), neighbors_NA, surface_voxels_NA, surface_neighbors_NA));
	}
	else {
	  cells.push_back(new image_cell(idx, meshes[idx], conf, comps.at(idx), 
					 neighbors[idx], cell_surface_voxels[idx], cell_surface_neighbors[idx]));
	}
	if(traj_ids.size() > 0) { 
	  cells.back()->trajectoryid = traj_ids.at(idx); 
	  if(traj_ids.at(idx) < 0) cells.back()->displayed = false; 
	  else cells.back()->displayed = true; 
	}
      }
    }

    // Build BVH on cells themselves. Intersect this BVH first before intersecting triangle mesh BVH.
    if(is_nuc) bvh_nuc = new BVH2<image_cell>(nuclei);
    else bvh = new BVH2<image_cell>(cells);
  }

  void image_cell::measure_dorsalfolds(vector<volume8*> &vs, vector<volume8*> &dispvols) {
    cout << "measure_dorsalfolds() trajectoryid = " << trajectoryid << endl;
    if(vs.size() < 5) return; // Make sure there are the right numebr of channels.
    float alpha = conf.voxel_alpha();

    dorsalfolds.total_surface_area = compute_sa();     // Total surface area

    //const uint8 bazooka_thresh = 75;
    float bazooka_thresh_percent = conf.baz_thresh_pcent;
    const int bazooka_channel = 5;  // smoothed Bazooka is in channel 5
    const int bazooka_channel_non_smooth = 2;
    const int par_channel = 1; // Par is always in channel 1.

    // Accumulate a histogram of bazooka values.
    vector<int> bazooka_hist(256, 0);
    volume8 &bazooka_smooth = *vs[bazooka_channel-1];
    for(int z = 0; z < bazooka_smooth.depth; z++)
      for(int y = 0; y < bazooka_smooth.height; y++)
	for(int x = 0; x < bazooka_smooth.width; x++) bazooka_hist[bazooka_smooth(x,y,z)]++;

    // The the threshold based on a percentage.
    int N_hist = 0; for(int i = 0; i < 256; i++) N_hist += bazooka_hist[i];
    int cum = 0; int bazooka_thresh = 0;
    while(bazooka_thresh < 256) {
      cum += bazooka_hist[bazooka_thresh];
      if(float(cum) / float(N_hist) >= bazooka_thresh_percent) break;
      bazooka_thresh++;
    }
						   
    mesh->sample(v_samples, *vs[bazooka_channel-1]);  // Sample median filtered volume along Bazooka labeling
    vector<int> above_thresh(v_samples.size(), 0);
    
    for(size_t vidx = 0; vidx < v_samples.size(); vidx++) 
      above_thresh[vidx] = v_samples[vidx] > bazooka_thresh ? 1 : 0;

    labels.clear();
    labels.push_back(above_thresh);     // Save labels for display.

    // Get centroids of triangles with intensity above threshold intensity.
    vector<vec3> troids;
    vector<triangle> tri;
    mesh->as_triangles(tri);
    for(size_t vidx = 0; vidx < above_thresh.size(); vidx++) {
      if(above_thresh[vidx] > 0) 
	troids.push_back(tri[vidx].centroid());
    }
    if(troids.size() == 0) { displayed = false; return; }
    // Compute average of those centroids to get center of plane.
    vec3 bazooka_mu(0,0,0);
    float N = troids.size();
    for(size_t t = 0; t < troids.size(); t++) bazooka_mu += troids[t] / N;

    // Compute PCA and determine which triangles are above the plane.
    vector<vec3> pca3; vector<float> pca3_x;
    PCA3(pca3, pca3_x, troids);

    // Get triangles above plane and below.
    plane3 plane; plane.from_pnorm(bazooka_mu, pca3[0]);
    vector<int> above_plane(v_samples.size(), 0);
    float pos_area = 0, neg_area = 0;
    for(size_t vidx = 0; vidx < tri.size(); vidx++) {
      triangle t = tri[vidx];
      t[0] *= alpha; t[1] *= alpha; t[2] *= alpha;
      if(plane.pos_side(tri[vidx].centroid())) { pos_area += t.area(); above_plane[vidx] = 1; }
      else { neg_area += t.area(); above_plane[vidx] = 2; }
    }
    dorsalfolds.above_baz_surface_area = min(pos_area, neg_area);
    labels.push_back(above_plane);

    volume8 &par = *vs[par_channel - 1]; 

    // Get voxel points in the mesh.
    vector<ivec3> vpts; mesh->volpts(vpts);
    double vxd = alpha, vxd3 = vxd * vxd * vxd;
    dorsalfolds.total_volume = vxd3 * (double)vpts.size();
    float pos_vol = 0, neg_vol = 0;
    float pos_par = 0, neg_par = 0; // Collect also par intensity above and below plane.
    for(size_t v = 0; v < vpts.size(); v++) {
      if(plane.pos_side(vpts[v].to_vec3())) {
	pos_vol += vxd3;
	pos_par += par(vpts[v]);
      }
      else {
	neg_vol += vxd3;
	neg_par += par(vpts[v]);
      }
    }
    // Use minimum volume.
    dorsalfolds.above_baz_volume = min(pos_vol, neg_vol);

    // Use volume to define apical and basal.
    bool pos_is_apical = pos_vol < neg_vol;
    if(pos_is_apical) dorsalfolds.basal_par_intensity = neg_par;
    else dorsalfolds.basal_par_intensity = pos_par;

    // Get total bazooka and par intensity with the cell.
    volume8 &bazooka = *vs[bazooka_channel_non_smooth - 1];
    dorsalfolds.total_par_intensity = 0;
    // Bazooka quantities.
    dorsalfolds.bazooka_nonpatch_intensity = 0;
    dorsalfolds.bazooka_nonpatch_volume = 0;
    dorsalfolds.bazooka_patch_intensity = 0;
    dorsalfolds.bazooka_patch_volume = 0;
    for(size_t v = 0; v < vpts.size(); v++) {
      if(vpts[v].x() < 0 || vpts[v].y() < 0 || vpts[v].z() < 0) continue;
      if(vpts[v].x() >= par.width || vpts[v].y() >= par.height || vpts[v].z() >= par.depth) continue;
      dorsalfolds.total_par_intensity += par(vpts[v]);
      if(bazooka(vpts[v]) > bazooka_thresh) {
	dorsalfolds.bazooka_patch_intensity += bazooka(vpts[v]);
	dorsalfolds.bazooka_patch_volume += vxd3;
      }
      else {
	dorsalfolds.bazooka_nonpatch_intensity += bazooka(vpts[v]);
	dorsalfolds.bazooka_nonpatch_volume += vxd3;
      }
    }
    vec3 cent = alpha * mesh->centroid;
    dorsalfolds.cent_x = cent[0]; dorsalfolds.cent_y = cent[1]; dorsalfolds.cent_z = cent[2];

    // Project vertices along the long dimension.
    vector<vec3> vec3vert(mesh->num_vert());
    for(int vidx = 0; vidx < mesh->num_vert(); vidx++) vec3vert[vidx] = alpha * (mesh->vertex(vidx) - mesh->centroid);
    pca3 = mesh->pca3;
    vec3 dir1 = pca3[2].max_project(vec3vert); // Make sure to rescale to micron axes.

    // Negate the long vector to get the other half.
    vec3 pca32_neg = -pca3[2];
    vec3 dir2 = pca32_neg.max_project(vec3vert);
    dorsalfolds.total_length = dir1.length() + dir2.length();

    vec3 pos1 = dir1 + alpha * mesh->centroid;
    vec3 pos2 = dir2 + alpha * mesh->centroid;

    // Larger/smaller z used to define apical.
    vec3 apical, basal;
    if(pos_is_apical) {
      if(plane.pos_side(pos1 / alpha)) { apical = pos1; basal = pos2; }
      else { apical = pos2; basal = pos1; }
    }
    else {
      if(plane.neg_side(pos1 / alpha)) { apical = pos1; basal = pos2; }
      else { apical = pos2; basal = pos1; }
    }
    //if(pos1[2] < pos2[2]) { apical = pos1; basal = pos2; } else { apical = pos2; basal = pos1; }

    vec3 basal_dir = basal - apical; basal_dir.normalize();
    vec3 bazooka_vec = basal_dir.proj(alpha * bazooka_mu - apical);
    dorsalfolds.apical_length = bazooka_vec.length();
    dorsalfolds.bazooka_pos = (dorsalfolds.apical_length /  dorsalfolds.total_length) * 100.0;

    cout << "par_dist 2 voxels in (microns)=" << 2 * alpha << endl;
    cout << "par_dist 4 voxels in (microns)=" << 4 * alpha << endl;

    vector<ivec3> vpts_par_4voxel; mesh->volpts(vpts_par_4voxel, 4);
    dorsalfolds.basal_4voxel_par_intensity = 0;
    dorsalfolds.basal_4voxel_par_volume = 0;
    for(size_t v = 0; v < vpts_par_4voxel.size(); v++) {
      if( (pos_is_apical  && plane.neg_side(vpts_par_4voxel[v].to_vec3())) ||
	  (!pos_is_apical && plane.pos_side(vpts_par_4voxel[v].to_vec3())) ) {
	dorsalfolds.basal_4voxel_par_intensity += par(vpts_par_4voxel[v]);
	dorsalfolds.basal_4voxel_par_volume += vxd3;
	//dispvols[0]->v(vpts_par_4voxel[v]) = 255;
      }
    }
    vector<ivec3> vpts_par_2voxel; mesh->volpts(vpts_par_2voxel, 2);
    dorsalfolds.basal_2voxel_par_intensity = 0;
    dorsalfolds.basal_2voxel_par_volume = 0;
    for(size_t v = 0; v < vpts_par_2voxel.size(); v++) {
      if( (pos_is_apical  && plane.neg_side(vpts_par_2voxel[v].to_vec3())) ||
	  (!pos_is_apical && plane.pos_side(vpts_par_2voxel[v].to_vec3())) ) {
	dorsalfolds.basal_2voxel_par_intensity += par(vpts_par_2voxel[v]);
	dorsalfolds.basal_2voxel_par_volume += vxd3;
	//dispvols[0]->v(vpts_par_2voxel[v]) = 190;
      }
    }
    mesh->boundary_edges(labels[label_idx], v0, v1);
  }

  void image_cell::measure_generic(vector<volume8*> &dispvols) {
    assign_neighbor_labels();
    surface_area = compute_sa();
    model_volume = mesh->compute_volume(conf.voxel_alpha());
    mesh->pca_dims(pca_dims, conf.voxel_alpha() );
    measure_lateral_bend();
    measure_neighbor_contact_properties();
    mesh->pca_longdir(apical_pos, basal_pos, 1.0);

    if(nucleus != NULL) {
      vec3 &nuc_centroid = nucleus->mesh->centroid;
      // NOTE: Using PCA axis from nucleus, may want to use cell mesh.
      vec3 nuc_longdir = nucleus->mesh->pca3.at(2);
      // Flip if pointing downward.
      if(nuc_longdir.z() > 0) nuc_longdir = -nuc_longdir;
      plane3 plane; plane.from_pnorm(nuc_centroid, nuc_longdir);
      vector<ivec3> above_pts, below_pts;
      mesh->volpts(above_pts, below_pts, plane);

      /*if(dispvols.size() > 0) {
	volume8 &dv = *dispvols.back();
	for(size_t a = 0; a < above_pts.size(); a++) 
	  dv(above_pts[a]) = 168;
	  }*/

      double vxd = conf.voxel_alpha(), vxd3 = vxd * vxd * vxd;
      vol_above_nuc = vxd3 * (double)above_pts.size();

      // Get surface area above nuc.
      sa_above_nuc = 0;
      vector<triangle> tri;
      mesh->as_triangles(tri);
      for(size_t tidx = 0; tidx < tri.size(); tidx++) {
	if( plane.pos_side(tri[tidx].centroid()) ) {
	  triangle t = tri[tidx];
	  t[0] *= vxd; t[1] *= vxd; t[2] *= vxd;
	  sa_above_nuc += t.area();
	}
      }

      // Get distance of nuc entroid from top.
      nuc_apical_dist = geom::distance3(vxd * apical_pos, vxd * nuc_centroid);
    }
  }

  struct fidxtroid_t {
    int fidx; vec3 troid; float dist;
    bool operator < (const fidxtroid_t &b) const { return dist < b.dist; }
  };
  struct sa_bend_t {
    float surface_area, bend_measure;
    bool operator < (const sa_bend_t &b) const { return surface_area > b.surface_area;  }
  };

  void image_cell::measure_lateral_bend() {
    if(mesh_face_neighbors.size() == 0) return;
    // Get unique labels in label vector.
    set<int> ulabels;
    for(size_t i = 0; i < mesh_face_neighbors.size(); i++) {
      if(mesh_face_neighbors[i] >= 0) ulabels.insert(mesh_face_neighbors[i]);    
    }
    vector<int> subidx;
    vector<vec3> vtable; vector<face> ftable;

    vector<sa_bend_t> sa_bend_vals;
    for(set<int>::iterator it = ulabels.begin(); it != ulabels.end(); it++) {
      // Get face indices of triangles with current label.
      subidx.clear();
      for(int fidx = 0; fidx < (int)mesh_face_neighbors.size(); fidx++) { if(mesh_face_neighbors[fidx] == *it) subidx.push_back(fidx); }
      // Build a mesh.
      mesh->as_mesh(subidx, vtable, ftable);
      geom::mesh *submesh = new geom::mesh(vtable, ftable);

      sa_bend_t sa_bend;
      sa_bend.bend_measure = submesh->bend_measure();
      sa_bend.surface_area = submesh->compute_sa(conf.voxel_alpha());
      sa_bend_vals.push_back(sa_bend);

      delete submesh;
    }    
    sort(sa_bend_vals.begin(), sa_bend_vals.end());

    max_cell_bend = 0;
    for(int s = 0; s < min(4, (int)sa_bend_vals.size()); s++) {
      max_cell_bend = max(max_cell_bend, sa_bend_vals[s].bend_measure);
    }
  }

  void image_cell::measure_neighbor_contact_properties() {
    neighbor_bend.clear();
    neighbor_contact_sa.clear();
    // Scan through neighbors.
    for(size_t n = 0; n < neighbors.size(); n++) {
      // Get faces in contact with the neighbor.
      vector<int> subidx;
      vector<vec3> vtable; vector<face> ftable;
      for(int fidx = 0; fidx < (int)mesh_face_neighbors.size(); fidx++) { 
	if(mesh_face_neighbors[fidx] == neighbors[n]) subidx.push_back(fidx);
      }
      // Build a mesh.
      mesh->as_mesh(subidx, vtable, ftable);
      geom::mesh *submesh = new geom::mesh(vtable, ftable);
      // Get contact surface area and bending.
      neighbor_bend.push_back(submesh->bend_measure());      
      neighbor_contact_sa.push_back(submesh->compute_sa(conf.voxel_alpha()));
      delete submesh;
    }
  }


  void image_stack::save_generic(string filename) {
    ofstream output(filename.c_str());
    output << "id\tsurface.area\tvolume" << endl;
    for(size_t i = 0; i < cells.size(); i++) {
      if(cells[i]->displayed) {
	output << cells[i]->idx + 1 << "\t" << cells[i]->surface_area << "\t" << cells[i]->model_volume << endl;
      }
    }
  }

  void image_cell::assign_neighbor_labels() {
    // Precompute the largest border needed since mesh might expand past surface voxels.

    // Get AABB for surface and mesh.
    vec3 surf_min, surf_max, mesh_min, mesh_max;
    geom::bounding_box(surf_min, surf_max, surface_voxels);
    mesh->bounding_box(mesh_min, mesh_max);

    // Find regions where mesh extends outside of surface.
    vec3 diff_min = surf_min - mesh_min, diff_max = mesh_max - surf_max;

    // Add a border to include those triangles.
    float brdr_diff = 0;
    for(int i = 0; i < 3; i++) {
      if(diff_min[i] > 0) brdr_diff = max(brdr_diff, diff_min[i]);
      if(diff_max[i] > 0) brdr_diff = max(brdr_diff, diff_max[i]);
    }
    //cout << "brdr_diff=" << ceil(brdr_diff) << endl;

    // fg = 255 at surface voxels
    volume8 surface(surface_voxels, 2 + ceil(brdr_diff));

    // Compute map to nearest surface voxel.
    volumeT<ivec3> nearest(surface.width, surface.height, surface.depth);
    surface.ComputeFT(nearest, 255);

    // Compute map from surface voxel to a neighbor.
    volume32 surf2neighbor(surface.width, surface.height, surface.depth);
    surf2neighbor.fill(-1);

    ivec3 shift_ivec(surface.x0, surface.y0, surface.z0);
    // NOTE: surf2nighbor has the labels 1-indexed.
    for(size_t v = 0; v < surface_neighbors.size(); v++) 
      surf2neighbor(surface_voxels[v] - shift_ivec) = surface_neighbors[v] + 1;

    // Shift triangles so they are within the surface volume.
    vec3 tri_shift(surface.x0, surface.y0, surface.z0);
    vector<triangle> tri; mesh->as_triangles(tri);
    for(size_t t = 0; t < tri.size(); t++) {
      for(int i = 0; i < 3; i++) tri[t][i] -= tri_shift;
    }

    vector<int> tri_labels(tri.size(), 0);     // Labels assigned to triangles.

    // Intersect triangles.
    for(size_t t = 0; t < tri.size(); t++) {
      vector<ivec3> res;
      nearest.tri_intersect(res, tri[t]);
      // Map intersected voxels to nearest surface voxel.
      for(size_t r = 0; r < res.size(); r++) {
	if(res[r][0] < 0 || res[r][0] >= nearest.width ||
	   res[r][1] < 0 || res[r][1] >= nearest.height ||
	   res[r][2] < 0 || res[r][2] >= nearest.depth) {
	  cerr << "before nearest error" << endl; exit(1);
	}
	res[r] = nearest(res[r]);
      }

      // Vote for a neighbor index.
      map<int, int> votes;
      for(size_t r = 0; r < res.size(); r++) {
	if(res[r][0] < 0 || res[r][0] >= surf2neighbor.width ||
	   res[r][1] < 0 || res[r][1] >= surf2neighbor.height ||
	   res[r][2] < 0 || res[r][2] >= surf2neighbor.depth) {
	  cerr << "after!!! nearest error" << endl; exit(1);
	}
	int neighbor_idx = surf2neighbor(res[r]);
	if(votes.find(neighbor_idx) == votes.end()) votes[neighbor_idx] = 1; else votes[neighbor_idx]++;
      }
      
      int max_idx = -1, max_cnt = 0;
      for(map<int, int>::iterator it = votes.begin(); it != votes.end(); it++) {
	int neighbor_idx = it->first;
	int neighbor_count = it->second;
	if(neighbor_count > max_cnt ) {
	  max_idx = neighbor_idx;
	  max_cnt = neighbor_count;
	}
      }
      tri_labels[t] = max_idx; // NOTE 1-indexed neighbors (subtract 1 to get index in cell_voxels!!!).
    }

    // NOTE: Use a planar surface to "cut" off cap regions.
    //vector<int> upd_tri_labels;
    //mesh->fill_bg(upd_tri_labels, tri_labels, 2, 0);
    //labels.push_back(upd_tri_labels);
    labels.clear();
    labels.push_back(tri_labels);

    // Subtract 1 so face corresponds to neighbor in cell_voxels. 
    mesh_face_neighbors = tri_labels;
    for(size_t m = 0; m < mesh_face_neighbors.size(); m++) mesh_face_neighbors[m] -= 1;
  }

  void image_stack::save_raw_segmentation(string filename, bool clean_shapes) {
    volume32 labels(iwidth, iheight, idepth);
    if(clean_shapes) {
      // TODO: Clean these shapes by dialation and erosion?
      // TODO: Also voxelize smoothed shapes.
      vector< vector<ivec3> > cleaned(cell_voxels.size()); 
      for(size_t c = 0; c < cell_voxels.size(); c++) {
	cout << "cleaning shape " << c << endl;
	volume8 v(cell_voxels[c], 2);
	v.fill_holes();
	v.dilate(1, 255, 0); v.dilate(1, 0, 255); // dilate and erode
	vector< vector<ivec3> > filled;
	v.components2(filled);
	if(filled.size() > 0) cleaned[c] = filled[0];
      }
      labels.seed(cleaned);
    }
    else labels.seed(cell_voxels); 

    labels.save_raw(filename);
  }

  void static_analysis_internal(vector<int> &seg_status, 
				vector<int> &assigned_truth,
				int &missed_truth,
				ivec3 dims, 
				float truth_alpha, 
				float seg_alpha,
				vector< vector<ivec3> > &truth_voxels, vector< vector<ivec3> > &seg_voxels) {
    // Label both segmented and truth regions.
    int w = dims[0], h = dims[1], d = dims[2];
    volume32 truth(w,h,d); truth.seed(truth_voxels);
    volume32 seg(w,h,d);   seg.seed(seg_voxels);

    vector<bool> split_truth(truth_voxels.size(), false);

    // Look for missed truth regions and truth regions that are split.
    missed_truth = 0;
    for(size_t t = 0; t < truth_voxels.size(); t++) {
      if(truth_voxels[t].size() == 0) continue;
      // Intersect and count overlap with segemnted regions.
      map<int,int> seg_idxs0; // seg region label -> number of overlap voxels
      vector<ivec3> &truth_vox = truth_voxels[t];
      for(size_t i = 0; i < truth_vox.size(); i++) {
	// Get overlaps for a truth region.
	if(seg(truth_vox[i]) > 0) {
	  if(seg_idxs0.find(seg(truth_vox[i])) == seg_idxs0.end()) seg_idxs0[seg(truth_vox[i])] = 1;
	  else seg_idxs0[seg(truth_vox[i])]++; 
	}
      } 
      // Make sure overlap is more than seg_alpha% of segmented region size and
      // more truth_alpha% of truth region size.
      int ok_overlap_cnt = 0;
      for(map<int,int>::iterator it = seg_idxs0.begin(); it != seg_idxs0.end(); it++) {
	// If overlap is > 10% of the nuc volume, save index.
	int seg_idx = it->first - 1; // -1 since labels are 1 indexed.
	int overlap = it->second;
	// Is the mutual overlap > 10% for each? 
	if( overlap > (int)(seg_voxels[seg_idx].size() * seg_alpha) && 
	    overlap > (int)(truth_vox.size() * truth_alpha) ) {
	  ok_overlap_cnt++;
	}
      }
      cout << "ok_overlap_cnt=" << ok_overlap_cnt << endl;
      // Insufficient overlap found, this is a missed region.
      if(ok_overlap_cnt == 0) missed_truth++;
      if(ok_overlap_cnt >= 2) split_truth[t] = true; // 2 or more regions, split truth region
    }

    // Detected: Single indexes found in region.
    // Merged: Multiple indexes found in region.
    // Junk: No  index found in region.
    // SplitMerge: Overlap with a split truth region and merged with another.
    // Split: Single  index, but that spans two or more regions.
    seg_status.resize(seg_voxels.size());
    assigned_truth.resize(seg_voxels.size());

    for(size_t s = 0; s < seg_voxels.size(); s++) {     // Scan through segmented regions.
      vector<ivec3> &seg_vox = seg_voxels[s];
      if(seg_vox.size() == 0) { 
	assigned_truth[s] = -1; seg_status[s] = analysis::StaticJunk; continue;
      }
      // Intersect and count overlap with truth regions.
      map<int,int> truth_idxs0; // truth region index -> number of overlap voxels
      for(size_t i = 0; i < seg_vox.size(); i++) {
	// Get overlap with truth regions.
	if(truth(seg_vox[i]) > 0) {
	  if(truth_idxs0.find(truth(seg_vox[i])) == truth_idxs0.end())  truth_idxs0[truth(seg_vox[i])] = 1;
	  else truth_idxs0[truth(seg_vox[i])]++; 
	}
      } 
      int ok_overlap_cnt = 0, split_overlap_cnt = 0;
      int putative_truth_idx = -1;
      for(map<int,int>::iterator it = truth_idxs0.begin(); it != truth_idxs0.end(); it++) {
	int truth_idx = it->first - 1;
	int overlap = it->second;
	// Is the mutual overlap OK? 
	//cout << overlap << "/" << seg_vox.size() * seg_alpha << "/" << truth_voxels[truth_idx].size() * truth_alpha  << endl;
	if(overlap > (int)(seg_vox.size() * seg_alpha) &&  overlap > (int)(truth_voxels[truth_idx].size() * truth_alpha) ) {
	  putative_truth_idx = truth_idx;
	  if(split_truth[truth_idx]) split_overlap_cnt++;  
	  else ok_overlap_cnt++;
	}
      }
      assigned_truth[s] = -1;
      if(ok_overlap_cnt == 1 && split_overlap_cnt == 0) {
	assigned_truth[s] = putative_truth_idx;
	seg_status[s] = analysis::StaticNoError; // Detected
      }
      else if(ok_overlap_cnt == 0 && split_overlap_cnt == 0) seg_status[s] = analysis::StaticJunk; // Junk.
      else if(ok_overlap_cnt == 0 && split_overlap_cnt == 1) seg_status[s] = analysis::StaticSplit; // Split
      else {
	if(split_overlap_cnt == 0) seg_status[s] = analysis::StaticMerge; // Merge
	else seg_status[s] = analysis::StaticSplitMerge;  // Merge and Split region.
      }
    }
  }

  // Determine which cells are poorly segmented using the nuclei channel for
  // each static time point. Analysis does not consider tracking.
  void image_stack::static_analysis_header(ofstream &output) {
    output << "nucs\t" << "cells\t" << "missed_nuc\t" << "detected\t" << "split\t" << "merged\t" << "junk\t" << "mergesplit";
  }
  void image_stack::static_analysis(ofstream &output) {
    // Remove cells and voxels that were truncated or touch the border region.
    ivec3 dims(iwidth, iheight, idepth);
    vector<int> seg_status, assigned_truth;
    int missed = 0;
    // >10% overlap of a truth region,  > 0.1% overlap of segmented region
    static_analysis_internal(seg_status, assigned_truth, missed, dims, 0.1, 0.001, nuc_voxels, cell_voxels);

    int detected_cell = 0, merge_cell = 0, junk_cell = 0, split_cell = 0, mergesplit_cell = 0;
    for(size_t s = 0; s < seg_status.size(); s++) {
      switch(seg_status[s]) {
      case 0: detected_cell++; break;
      case 1: 
	if(is_border_touch(cell_voxels.at(s)) == false) junk_cell++; 
	break;
      case 2: split_cell++; break;
      case 3: merge_cell++; break;
      case 4: mergesplit_cell++; break;
      }
    }

    output << nuc_voxels.size() << '\t'; output << cell_voxels.size() << '\t';
    output << missed << '\t';            output << detected_cell << '\t';
    output << split_cell << '\t';        output << merge_cell << '\t';
    output << junk_cell << '\t';         output << mergesplit_cell ;
  }

  void associate_motion3(vector<vec3> &cur_cents, vector<vec3> &cur_motion, 
			 vector<int> &cur_match, 
			 vector<vec3> &nxt_cents, float nearestk_dist);

  // Centroid tracking with a motion model.
  void hyperstack::track_nucs_centroid2() {
    cout << "track_nucs_centroid2()" << endl;
    if(stacks.size() == 0) return;

    // Assign initial trajectory IDs.
    vector<int> &init_traj_ids = stacks[0]->nuc_traj_ids;  
    init_traj_ids.resize(stacks[0]->nuc_voxels.size());
    stacks[0]->nuc_motion.resize(stacks[0]->nuc_voxels.size());
    for(size_t t = 0; t < init_traj_ids.size(); t++) {
      init_traj_ids[t] = t;
      stacks[0]->nuc_motion[t] = vec3(0,0,0);
    }

    unsigned int last_traj_id = init_traj_ids.size(); // Save the last used trajectory ID.

    for(size_t t = 0; t < stacks.size() - 1; t++) {
      image_stack *cur_stack = stacks[t], *nxt_stack = stacks[t+1];   // Get current and next frme.

      vector<vec3> nxt_nucs; nxt_stack->nuc_centroids(nxt_nucs);      
      vector<vec3> cur_nucs; cur_stack->nuc_centroids(cur_nucs);
      vector<vec3> &cur_motion = cur_stack->nuc_motion;

      // Build map from cur_idx to assigne nxt_idx in nxt_nucs.
      vector<int> cur_match(cur_nucs.size(), -1); 

      vector<vec3> no_motion(cur_nucs.size());
      for(size_t i = 0; i < no_motion.size(); i++) no_motion[i] = vec3(0,0,0);
      associate_motion3(cur_nucs, no_motion, cur_match, nxt_nucs, conf.nearestk_dist);

      vector<int> cur_match_motion(cur_nucs.size(), -1);
      associate_motion3(cur_nucs, cur_motion, cur_match_motion, nxt_nucs, conf.nearestk_dist);

      int matches_killed = 0;
      for(int cur_idx = 0; cur_idx < (int)cur_match.size(); cur_idx++) {
	// Disagreeing assignment between motion and nearest.
	if(cur_match[cur_idx] != cur_match_motion[cur_idx]) {
	  cur_match[cur_idx] = -1; matches_killed++;
	}
      }
      cout << "matches_killed=" << matches_killed << endl;
      // Match on motion only.
      associate_motion3(cur_nucs, cur_motion, cur_match, nxt_nucs, conf.nearestk_dist);

      // Match on nearest.
      for(size_t c = 0; c < cur_motion.size(); c++) cur_motion[c] = vec3(0,0,0);
      associate_motion3(cur_nucs, cur_motion, cur_match, nxt_nucs, conf.nearestk_dist);

      // TODO: Rematch using motion only, not killing previous assignments.
      vector<vec3> &nxt_motion = nxt_stack->nuc_motion;
      nxt_motion.resize(nxt_nucs.size());
      for(size_t n = 0; n < nxt_motion.size(); n++) nxt_motion[n] = vec3(0,0,0);

      vector<int> &cur_traj_ids = cur_stack->nuc_traj_ids;
      vector<int> &nxt_traj_ids = nxt_stack->nuc_traj_ids;

      nxt_traj_ids.resize(nxt_nucs.size());
      for(int nxt_idx  = 0; nxt_idx < (int)nxt_traj_ids.size(); nxt_idx++) nxt_traj_ids[nxt_idx] = -1;

      int lost_track = 0;
      for(int cur_idx = 0; cur_idx < (int)cur_nucs.size(); cur_idx++) {
	if(cur_match[cur_idx] == -1)  { lost_track++; continue; }
	// Transfer trajectory ID.
	nxt_traj_ids.at(cur_match[cur_idx]) = cur_traj_ids[cur_idx];
	nxt_motion[cur_match[cur_idx]] = nxt_nucs.at(cur_match[cur_idx]) - cur_nucs[cur_idx];
      }
      cout << "lost_track=" << lost_track << endl;
      // Assign new trajectory IDs.
      int new_ids = 0;
      for(int nxt_idx = 0; nxt_idx < (int)nxt_nucs.size(); nxt_idx++) {
	if(nxt_traj_ids[nxt_idx] == -1) {
	  nxt_traj_ids[nxt_idx] = last_traj_id++;
	  new_ids++;
	}
      }
      cout << "new_ids=" << new_ids << endl;
    }    

    // Compute nuclei trajectory lengths.
    nuc_traj_lengths.resize(last_traj_id);
    for(size_t i = 0; i < nuc_traj_lengths.size(); i++) nuc_traj_lengths[i] = 0;

    for(size_t t = 0; t < stacks.size(); t++) {
      image_stack *cur_stack = stacks[t];
      vector<int> &nuc_traj_ids = cur_stack->nuc_traj_ids;
      map<int, int> &nuc_traj_id2nuc = cur_stack->nuc_traj_id2nuc;
      for(size_t i = 0; i < nuc_traj_ids.size(); i++) {
	nuc_traj_lengths.at(nuc_traj_ids[i])++;
	nuc_traj_id2nuc[nuc_traj_ids[i]] = i;
      }
    }

    cout << "track_nucs_centroid2() *DONE*" << endl;
  }

  // Nearest index and distance squared to the index.
  struct idxdist_t { 
    int next_idx; float distsq; 
    bool operator < (const idxdist_t &b) const { return distsq < b.distsq; }
  };

  // Get indicies of nearest K points to P.
  void NearestK(vector<idxdist_t> &nearestK, vec3 p, float dthresh, vector<vec3> &vs) {
    // Get index and distance.
    float dthreshsq = dthresh * dthresh;
    vector<idxdist_t> idxdist(vs.size());
    for(int idx = 0; idx < (int)vs.size(); idx++) {
      idxdist[idx].next_idx = idx;  
      idxdist[idx].distsq = geom::distance3sq(p, vs[idx]);
    }
    // Sort from closest to furthest and save those within a distance threshold
    sort(idxdist.begin(), idxdist.end());
    nearestK.clear();
    for(int i = 0; i < (int)idxdist.size(); i++) {
      if(idxdist[i].distsq <= dthreshsq) nearestK.push_back(idxdist[i]);
      else break;
    }
  }


  void associate_motion3(vector<vec3> &cur_cents, vector < vector<idxdist_t> > &cur_matches,
			 vector<int> &cur_match,  vector<vec3> &nxt_cents) {

    vector<int> nxt_assign(nxt_cents.size(), -1);      // Index in cur_cents that is assigned to nxt_nuc.
    vector<int> nxt_assign_prev(nxt_cents.size(), -1); // Index in cur_cents previously assigned to nxt_nuc.

    // Make sure number of centroids is the same as mapping from cell -> match in next.
    if(cur_match.size() != cur_cents.size()) { cerr << "cur_match.size() != cur_cents.size()" << endl; exit(1); } 
    if(cur_matches.size() != cur_cents.size()) { cerr << "cur_matches.size() != cur_cents.size()" << endl; exit(1); }

    // Get previous assignments.
    for(int cur_idx = 0; cur_idx < (int)cur_cents.size(); cur_idx++) {
      if(cur_match[cur_idx] != -1) {
	nxt_assign_prev.at(cur_match[cur_idx]) = cur_idx; // prevents improving on match to this centroid
	nxt_assign.at(cur_match[cur_idx]) = cur_idx;
      }
    }
    
    vector<int> cur_matchidx(cur_cents.size(), 0); // Index of the match to use in found matches, for a current cell.
    // Handle cases where no matches availible and previously assigned.
    for(int cur_idx = 0; cur_idx < (int)cur_cents.size(); cur_idx++) {
      if(cur_match[cur_idx] != -1) cur_matchidx[cur_idx] = -1;  // Skip and set match index to -1 so not processed below.
      else if(cur_matches[cur_idx].size() == 0) cur_matchidx[cur_idx] = -1;
    }

    // Scan through cells, greedily assigning on smallest distance. 
    int conflicts = cur_cents.size();
    while(conflicts > 0) {
      conflicts = 0;
      for(int cur_idx = 0; cur_idx < (int)cur_cents.size(); cur_idx++) {
	// Skip those with cells with matches availible or were already assigned.
	if(cur_matchidx[cur_idx] == -1) continue;

	// Get current match nxt index.
	int nxt_idx = cur_matches[cur_idx][cur_matchidx[cur_idx]].next_idx;

	if(nxt_assign[nxt_idx] == -1) nxt_assign[nxt_idx] = cur_idx; // Easy case, nothing assigned, assign match.
	else if(nxt_assign[nxt_idx] == cur_idx) continue; // No conflict.
	else {
	  // Conflict found.
	  conflicts++;
	  if(nxt_assign_prev[nxt_idx] != -1) { 	 // Conflict with nxt nuc index that was previously assigned.
	    cur_matchidx[cur_idx]++; // Try another nxt nuc index in topK set.	    
	    if(cur_matchidx[cur_idx] >= (int)cur_matches[cur_idx].size()) cur_matchidx[cur_idx] = -1;
	  }
	  else {
	    // Otherwise, try to improve assignment.
	    int old_idx = nxt_assign.at(nxt_idx);
	    float old_dist = cur_matches[old_idx][cur_matchidx[old_idx]].distsq;
	    float new_dist = cur_matches[cur_idx][cur_matchidx[cur_idx]].distsq;

	    if( old_dist < new_dist ) {
	      // The previous assignment is better, so keep and try next best match for current cell.
	      cur_matchidx[cur_idx]++;
	      if(cur_matchidx[cur_idx] >= (int)cur_matches[cur_idx].size()) cur_matchidx[cur_idx] = -1;
	    }
	    else {
	      // Match imrpoves the distance.
	      nxt_assign[nxt_idx] = cur_idx;// Use the current match instead.
	      cur_matchidx[old_idx]++;  // Try a better match for the old index.
	      if(cur_matchidx[old_idx] >= (int)cur_matches[old_idx].size()) cur_matchidx[old_idx] = -1;
	    }
	  }
	}
      }	
    } 

    for(int cur_idx = 0; cur_idx < (int)cur_cents.size(); cur_idx++) {
      if(cur_match[cur_idx] != -1) continue; // Already assigned.
      else {
	// Otherwise, use new nxt_nuc assignment, if it exists.
	if(cur_matchidx[cur_idx] != -1) {
	  cur_match[cur_idx] = cur_matches[cur_idx][cur_matchidx[cur_idx]].next_idx;
	}
      }
    }
  }

  void associate_motion3(vector<vec3> &cur_cents, vector<vec3> &cur_motion, 
			 vector<int> &cur_match, 
			 vector<vec3> &nxt_cents, float nearestk_dist) {    
    vector < vector<idxdist_t> > cur_matches_3d(cur_cents.size());
    // Get nearest K within a given distance. 
    for(int cur_idx = 0; cur_idx < (int)cur_cents.size(); cur_idx++) {
      NearestK(cur_matches_3d[cur_idx], cur_cents[cur_idx] + cur_motion[cur_idx], nearestk_dist, nxt_cents);
    }
    // Get 3d matches.
    associate_motion3(cur_cents, cur_matches_3d, cur_match, nxt_cents);
  }


  // Runs fractional Hausdorff matching in separate threads. 
  struct match_thread : public QThread {
    volume32 *curEDT;
    vector< vector<ivec3> > &nxt_onucs;
    vector<float> &nxt_dists;
    vector<ivec3> &nxt_shifts;
    int nuc_split_dist;
    int thread_start, thread_step;
    match_thread(volume32 *curEDT,
		 vector< vector<ivec3> > &nxt_onucs,
		 vector<float> &nxt_dists,
		 vector<ivec3> &nxt_shifts, 
		 int nuc_split_dist,
		 int thread_start, int thread_step) :
      curEDT(curEDT),
      nxt_onucs(nxt_onucs), nxt_dists(nxt_dists), nxt_shifts(nxt_shifts), 
      nuc_split_dist(nuc_split_dist),
      thread_start(thread_start), thread_step(thread_step) { }
    void run() {
      for(size_t n = thread_start; n < nxt_onucs.size(); n += thread_step)
	nxt_dists[n] = curEDT->frac_hausdorff_match(nxt_shifts[n], nxt_onucs[n], nuc_split_dist, 0.95);
    }
  };

  // Centroid tracking with a motion model.
  void hyperstack::track_cells_centroid2() {
    if(stacks.size() == 0) return;

    // Assign initial trajectory IDs.
    vector<int> &init_traj_ids = stacks[0]->traj_ids; 
    init_traj_ids.resize(stacks[0]->cell_voxels.size());

    // Set motion model vector to (0,0,0).
    stacks[0]->cell_motion.resize(stacks[0]->cell_voxels.size());
    for(size_t t = 0; t < init_traj_ids.size(); t++) {
      init_traj_ids[t] = t;
      stacks[0]->cell_motion[t] = vec3(0,0,0);
    }

    unsigned int last_traj_id = init_traj_ids.size(); // Save the last used trajectory ID.

    for(size_t t = 0; t < stacks.size() - 1; t++) {
      image_stack *cur_stack = stacks[t], *nxt_stack = stacks[t+1];   // Get current and next frme.

      // Get (EDT weighted) cell centroids and motion data.
      vector<vec3> nxt_cells; nxt_stack->cell_centroids(nxt_cells);      
      vector<vec3> cur_cells; cur_stack->cell_centroids(cur_cells);
      vector<vec3> &cur_motion = cur_stack->cell_motion;

      // Build map from cur_idx to assign nxt_idx in nxt_cells.
      vector<int> cur_match(cur_cells.size(), -1); 

      // Find nearest gready match, no motion.
      vector<vec3> no_motion(cur_cells.size());
      for(size_t i = 0; i < no_motion.size(); i++) no_motion[i] = vec3(0,0,0);
      associate_motion3(cur_cells, no_motion, cur_match, nxt_cells, conf.nearestk_dist);

      // Find nearest greedy match using motion model.
      vector<int> cur_match_motion(cur_cells.size(), -1);
      associate_motion3(cur_cells, cur_motion, cur_match_motion, nxt_cells, conf.nearestk_dist);

      // Kill cases of disagreemnt.
      for(int cur_idx = 0; cur_idx < (int)cur_match.size(); cur_idx++) {
	// Disagreeing assignment between motion and nearest.
	if(cur_match[cur_idx] != cur_match_motion[cur_idx]) cur_match[cur_idx] = -1;
      }
      // Rematch on motion only.
      associate_motion3(cur_cells, cur_motion, cur_match, nxt_cells, conf.nearestk_dist);

      // Rematch using no motion at all.
      associate_motion3(cur_cells, no_motion, cur_match, nxt_cells, conf.nearestk_dist);

      vector<vec3> &nxt_motion = nxt_stack->cell_motion;
      nxt_motion.resize(nxt_cells.size());
      for(size_t n = 0; n < nxt_motion.size(); n++) nxt_motion[n] = vec3(0,0,0);

      vector<int> &cur_traj_ids = cur_stack->traj_ids;
      vector<int> &nxt_traj_ids = nxt_stack->traj_ids;

      nxt_traj_ids.resize(nxt_cells.size());
      for(int nxt_idx  = 0; nxt_idx < (int)nxt_traj_ids.size(); nxt_idx++) nxt_traj_ids[nxt_idx] = -1;

      int lost_track = 0, num_tracked = 0;
      for(int cur_idx = 0; cur_idx < (int)cur_cells.size(); cur_idx++) {
	if(cur_match[cur_idx] == -1)  { lost_track++; continue; }
	// Transfer trajectory ID.
	nxt_traj_ids.at(cur_match[cur_idx]) = cur_traj_ids[cur_idx];
	// Update motion model.
	nxt_motion[cur_match[cur_idx]] = 0.5 * (nxt_cells.at(cur_match[cur_idx]) - cur_cells[cur_idx]) + 0.5 * cur_motion[cur_idx];
	num_tracked++;
      }
      cout << "lost_track=" << lost_track << endl;
      cout << "num_tracked=" << num_tracked << endl;
      // Assign new trajectory IDs.
      int new_ids = 0;
      for(int nxt_idx = 0; nxt_idx < (int)nxt_cells.size(); nxt_idx++) {
	if(nxt_traj_ids[nxt_idx] == -1) {
	  nxt_traj_ids[nxt_idx] = last_traj_id++;
	  new_ids++;
	}
      }
      cout << "new_ids=" << new_ids << endl;
    }    
    // ------------------------------------------ tracking is done

    // Get trajectory lengths.
    traj_lengths.resize(last_traj_id);
    for(size_t i = 0; i < traj_lengths.size(); i++) traj_lengths[i] = 0;
    for(size_t t = 0; t < stacks.size(); t++) {
      image_stack *cur_stack = stacks[t];
      vector<int> &traj_ids = cur_stack->traj_ids;
      // Save mapping from trajectory ID to cell shape ID.
      map<int, int> &traj_id2cell = cur_stack->traj_id2cell;
      for(size_t i = 0; i < traj_ids.size(); i++) {
	if(traj_ids[i] >= 0) {
	  traj_lengths.at(traj_ids[i])++;
	  traj_id2cell[traj_ids[i]] = i;
	}
      }
    }
  }

  // Assigns nuc to a cell shape based on overlap
  void image_stack::assign_nucs() {
    cout << "assign_nucs()" << endl;
    // Clear assigned nuclei.
    assigned_nucs.resize(cell_voxels.size());
    for(size_t a = 0; a < assigned_nucs.size(); a++) assigned_nucs[a] = -1;

    vector<int> assigned_overlap(assigned_nucs.size(), 0);

    volume32 cell_labels(iwidth, iheight, idepth); 
    cell_labels.seed(cell_voxels);

    for(size_t n = 0; n < nuc_voxels.size(); n++) {
      vector<ivec3> &vox = nuc_voxels[n];  
      // Get overlapping cells and keep track of overlap area.
      map<int, int> cell2olap;
      for(size_t i = 0; i < vox.size(); i++)  { 
	int cell_idx = cell_labels(vox[i]) - 1;
	if(cell_idx >= 0) {
	  if(cell2olap.find(cell_idx) != cell2olap.end()) cell2olap[cell_idx]++; 
	  else cell2olap[cell_idx] = 1;
	}
      }
      // Assign nuc to overlapping cell with greatest overlap
      int max_olap = 0, max_cell_idx = -1;
      for(map<int,int>::iterator it = cell2olap.begin(); it != cell2olap.end(); it++) {
	int cell_idx = it->first;
	int overlap = it->second;
	// Minimum overlap must be > 10% of the nuc size.
	if(overlap > (int)vox.size() / 10) {
	  if(overlap > max_olap) {
	    max_olap = overlap;
	    max_cell_idx = cell_idx;
	  }
	}
      }
      if(max_cell_idx >= 0) {
	// NOTE: A merged cell might be assigned two nuclei.
	if(assigned_nucs[max_cell_idx] == -1) {
	  assigned_nucs[max_cell_idx] = n;
	  assigned_overlap[max_cell_idx] = max_olap;
	}
	else if( assigned_overlap[max_cell_idx] < max_olap ) { 
	  cout << "assigned_overlap[max_cell_idx] < max_olap" << endl;
	  //exit(1);
	  assigned_nucs[max_cell_idx] = n;
	  assigned_overlap[max_cell_idx] = max_olap;
	}
      }
    }
  }

  // Slices voxels through an orthogonal slice of the cell segmentation at dim_pos.
  void image_stack::ortho_comps_cells(vector< vector<ivec3> > &comp, int dim, int dim_pos) {
    if(!(dim == 0 || dim == 1 || dim == 2)) { cerr << "ortho_comps():dim != 0,1,2" << endl; exit(1); }
    comp.resize(cell_voxels.size());
    for(size_t c = 0; c < cell_voxels.size(); c++) {
      vector<ivec3> &vox = cell_voxels[c];
      vector<ivec3> &slice = comp[c];
      slice.clear();
      for(size_t i = 0; i < vox.size(); i++) { if(vox[i][dim] == dim_pos) slice.push_back(vox[i]);  }
    }
  }
  
  // Repeats the same but for nuc voxels.
  void image_stack::ortho_comps_nucs(vector< vector<ivec3> > &comp, int dim, int dim_pos) {
    if(!(dim == 0 || dim == 1 || dim == 2)) { cerr << "ortho_comps():dim != 0,1,2" << endl; exit(1); }
    comp.resize(nuc_voxels.size());
    for(size_t c = 0; c < nuc_voxels.size(); c++) {
      vector<ivec3> &vox = nuc_voxels[c];
      vector<ivec3> &slice = comp[c];
      slice.clear();
      for(size_t i = 0; i < vox.size(); i++) { if(vox[i][dim] == dim_pos) slice.push_back(vox[i]);  }
    }
  }

  bool image_stack::is_border_touch(vector<ivec3> &vox) {
    bool border_touch = false;
    int brdr = clear_border() + 1; 
    for(size_t i = 0; i < vox.size(); i++) {
      if(vox[i][0] <= brdr || vox[i][1] <= brdr || vox[i][2] <= brdr ||
	 vox[i][0] >= iwidth - 1 - brdr || vox[i][1] >= iheight - 1 - brdr || vox[i][2] >= idepth - 1 - brdr) {
	border_touch = true; break;
      }
    }
    return border_touch;
  }


  struct seg_thread  : public QThread {
    bool split_nucs;
    int time_start, time_step;
    edgeconfig &conf;
    vector<image_stack*> &stacks;
    seg_thread(bool split_nucs, int time_start, int time_step, edgeconfig &conf, vector<image_stack *>  &stacks) : 
      split_nucs(split_nucs), time_start(time_start), time_step(time_step), conf(conf), stacks(stacks) {  }
    void run() {
      // Move this to a thread.
      for(size_t t = time_start; t < stacks.size(); t += time_step) {
	// Create a binary volume with edge voxels.
	volume8 bv(stacks[t]->iwidth, stacks[t]->iheight, stacks[t]->idepth);
	vector<ivec3> &edge_voxels = stacks[t]->edge_voxels;
	for(size_t e = 0; e < edge_voxels.size(); e++) bv(edge_voxels[e]) = 255;
	if(conf.analysis_id == analysis::NucsMembranes || conf.analysis_id == analysis::NucsOnly || split_nucs) {	    
	  // Use hybrid edge/nuc segmentatiion to generate cells.
	  stacks[t]->cell_voxels.clear();
	  HSegment3_Nucs(*stacks[t]->vs.at(conf.edge_channel-1), bv, stacks[t]->nuc_voxels,
			 1, conf.max_hole_rad, conf.internal_limit, 
			 stacks[t]->cell_voxels, 
			 conf.min_comp_vol, conf.max_comp_vol, conf.noise_comp_vol);
	  if(split_nucs) {
	    // Use segemntatiion to split cells.
	    stacks[t]->split_nucs();
	    stacks[t]->cell_voxels.clear();
	  }
	}
	else {
	  vector< vector<ivec3> > nuc_vox_empty;
	  HSegment3_Nucs(*stacks[t]->vs.at(conf.edge_channel-1), bv, nuc_vox_empty,
			 1, conf.max_hole_rad, conf.internal_limit, 
			 stacks[t]->cell_voxels, 
			 conf.min_comp_vol, conf.max_comp_vol, conf.noise_comp_vol);	     
	}
	stacks[t]->find_cell_neighbors2();
       }
     }
   };

  void hyperstack::run_segmentation(bool split_nucs) {
    vector<seg_thread *> threads(conf.hs_threads);
    for(int t = 0; t < conf.hs_threads; t++) threads[t] = new seg_thread(split_nucs, t, conf.hs_threads, conf, stacks);
    for(int t = 0; t < conf.hs_threads; t++) threads[t]->start();
    for(int t = 0; t < conf.hs_threads; t++) threads[t]->wait();
    for(int t = 0; t < conf.hs_threads; t++) delete threads[t];
  }

  // Performs single stack processing in a thread.
  struct proc_thread  : public QThread {
    int time_start, time_step;
    edgeconfig &conf;
    vector<image_stack*> &stacks;
    proc_thread(int time_start, int time_step, edgeconfig &conf, vector<image_stack *>  &stacks) : 
      time_start(time_start), time_step(time_step), conf(conf), stacks(stacks) { }
    void run() { for(size_t t = time_start; t < stacks.size(); t += time_step) stacks[t]->process(); }
  };

  // Builds cell shapes in separate threads.
  struct build_thread : public QThread {
    int time_start, time_step;
    vector<image_stack*> &stacks;
    build_thread(int time_start, int time_step, vector<image_stack *>  &stacks) : 
      time_start(time_start), time_step(time_step), stacks(stacks) { }
    void run() {
      for(size_t t = time_start; t < stacks.size(); t+=time_step) {
	// TODO: Interface to build_meshes is fucked up. Fix and clean up.
	stacks[t]->build_meshes(stacks[t]->cell_voxels, stacks[t]->cell_neighbors);
	vector< vector<int> > neighbors_NA;
	stacks[t]->build_meshes(stacks[t]->nuc_voxels, neighbors_NA, true);
	stacks[t]->transfer_nuc_meshes();
      }
    }
  };

  // Loads .TIFs as 16-bit images.
  void hyperstack::load_tiff_stacks() {
    // TODO: Need to handle out of range and errors better at load.
    // Loads each stack in the hyper stack.
    for(size_t fid = 0; fid < files.size(); fid++) {
      for(int frame = conf.startframe; frame <= conf.endframe; frame += conf.stepframe) {
	cout << "frame=" << frame << endl;
	cout << "conf.channels=" << conf.channels << endl;
	vector<volume16 *> vs16(conf.channels, NULL);
	cout << files[fid] << " " << conf.slices << " " << frame << " " << conf.channels << endl;
	int failed_load = 0;
	for(int c = 0; c < conf.channels; c++) vs16[c] = new volume16;
	for(int c = 0; c < conf.channels; c++) {
	  if(!vs16[c]->load(files[fid], conf.slices, conf.channels, frame, c)) failed_load++;
	}
	cout << "done loading" << endl;
	if(failed_load == 0) {
	  image_stack *stack = new image_stack(basename, vs16, conf);
	  cout << "creating stack" << endl;
	  stack->frame_num = frame;
	  if(files.size() > 1) stack->frame_num = fid;
	  stacks.push_back(stack);
	}
	else {
	  cout << "Failed to load " << failed_load << " channels." << endl;
	  for(int c = 0; c < conf.channels; c++) delete vs16[c];
	}
      }
    }
  }


  void hyperstack::process() {
    load_tiff_stacks();

    // Run single frame processing.
    vector<proc_thread *> threads(conf.hs_threads);
    for(int t = 0; t < conf.hs_threads; t++) threads[t] = new proc_thread(t, conf.hs_threads, conf, stacks);
    for(int t = 0; t < conf.hs_threads; t++) threads[t]->start();
    for(int t = 0; t < conf.hs_threads; t++) threads[t]->wait();
    for(int t = 0; t < conf.hs_threads; t++) delete threads[t];
    threads.clear();

    //if(conf.analysis_id == analysis::SPIMTest) {
    //stacks[0]->process_align(stacks[1]);
    //stacks[2]->process_align(stacks[3]);
    //}
    
    if(conf.analysis_id == analysis::NucsOnly || 
       conf.analysis_id == analysis::NucsMembranes || 
       conf.analysis_id == analysis::ManualGroundTruth ) {
      if(conf.repair_nucs) run_segmentation(true);
      track_nucs_centroid2();  
    }
    if(conf.analysis_id == analysis::StaticAnalysis) {
      if(conf.repair_nucs) run_segmentation(true);
    }
    if(conf.analysis_id == analysis::NucsOnly || conf.analysis_id == analysis::SPIMTest) {
      // Just build 3d models of nuclei.
      for(size_t t = 0; t < stacks.size(); t++) {
	stacks[t]->traj_ids = stacks[t]->nuc_traj_ids;
	stacks[t]->traj_id2cell = stacks[t]->nuc_traj_id2nuc;
	if(conf.run_mesh) {
	  vector< vector<int> > neighbors_NA;
	  stacks[t]->build_meshes(stacks[t]->nuc_voxels, neighbors_NA, false);
	}
      }
    }
    if(conf.analysis_id == analysis::NucsMembranes || conf.analysis_id == analysis::CellShapes     ||
       conf.analysis_id == analysis::DorsalFolds   || conf.analysis_id == analysis::StaticAnalysis ||
       conf.analysis_id == analysis::ManualGroundTruth) {
      if(conf.run_mesh) {
	run_segmentation();
	if(conf.analysis_id == analysis::NucsMembranes) {
	  // Transfer trajectory IDs from nuclei.
	  for(size_t t = 0; t < stacks.size(); t++) stacks[t]->transfer_traj_ids();

	  // Initialize trajectory lengths.
	  traj_lengths.resize(nuc_traj_lengths.size());
	  for(size_t i = 0; i < traj_lengths.size(); i++) traj_lengths[i] = 0;

	  // Compute trajectory lengths and keep mapping from trajectory 2 cell index.
	  for(size_t t = 0; t < stacks.size(); t++) {
	    image_stack *cur_stack = stacks[t];
	    vector<int> &traj_ids = cur_stack->traj_ids;
	    map<int, int> &traj_id2cell = cur_stack->traj_id2cell;
	    for(size_t i = 0; i < traj_ids.size(); i++) {
	      if(traj_ids[i] >= 0) {
		if(traj_id2cell.find(traj_ids[i]) != traj_id2cell.end() ) {
		  cerr << "double assignment of a trajectory" << endl;
		  exit(1);
		}
		traj_lengths.at(traj_ids[i])++;
		traj_id2cell[traj_ids[i]] = i;
	      }
	    }
	  }
	}
	else { 
	  cout << "running track cells centroid 2" << endl;
	  track_cells_centroid2();
	}
	cout << "building meshes" << endl;
	// Build meshes for cells and nucs. Threaded implemenation of the following:
	/*for(size_t t = 0; t < stacks.size(); t++) stacks[t]->build_meshes(stacks[t]->cell_voxels);
	  for(size_t t = 0; t < stacks.size(); t++) stacks[t]->build_meshes(stacks[t]->nuc_voxels, true);*/
	vector<build_thread *> threads2(conf.hs_threads);
	for(int t = 0; t < conf.hs_threads; t++) threads2[t] = new build_thread(t, conf.hs_threads, stacks);
	for(int t = 0; t < conf.hs_threads; t++) threads2[t]->start();
	for(int t = 0; t < conf.hs_threads; t++) threads2[t]->wait();
	for(int t = 0; t < conf.hs_threads; t++) delete threads2[t];
	threads2.clear();
      }
    }
    cout << "done hyperstack process" << endl;
  }


  struct measure_traj_thread : public QThread {
    edgeconfig &conf;
    vector<image_stack *> &stacks;
    vector<int> &traj_lengths;
    int traj_start, traj_step;
    measure_traj_thread(edgeconfig &conf, vector<image_stack *> &stacks, 
			vector<int> &traj_lengths, int traj_start, int traj_step) : 
      conf(conf),
      stacks(stacks),
      traj_lengths(traj_lengths),
      traj_start(traj_start), traj_step(traj_step)  { }
    void run() {
      for(int traj_id = traj_start; traj_id < (int)traj_lengths.size(); traj_id += traj_step) {
	if(traj_lengths[traj_id] <= 0) continue;
	cout << "traj_id=" << traj_id << endl;
	for(size_t t = 0; t < stacks.size(); t++) {
	  map<int,int> &traj_id2cell = stacks[t]->traj_id2cell;
	  if(traj_id2cell.find(traj_id) == traj_id2cell.end()) continue;
	  image_cell *cell = stacks[t]->cells.at(traj_id2cell[traj_id]);
	  if(conf.analysis_id == analysis::DorsalFolds)  cell->measure_dorsalfolds(stacks[t]->vs, stacks[t]->dispvols);
	  else cell->measure_generic(stacks[t]->dispvols);
	}	
      }
    }
  };

  void hyperstack::measure_trajectories_threaded() {
    // Apply measurements to cells using given number of CPU threads.
    vector<measure_traj_thread *> threads(conf.hs_threads);
    for(int t = 0; t < conf.hs_threads; t++) threads[t] = new measure_traj_thread(conf, stacks, traj_lengths, t, conf.hs_threads);
    for(int t = 0; t < conf.hs_threads; t++) threads[t]->start();
    for(int t = 0; t < conf.hs_threads; t++) threads[t]->wait();
    for(int t = 0; t < conf.hs_threads; t++) delete threads[t];
  }

  void hyperstack::save_contact_analysis(string filename) {
    measure_trajectories_threaded();

    ofstream output(filename.c_str());
    output << "frame\t"; // 1
    output << "time\t"; // 2 
    output << "trajectory.id\t";  // 3
    output << "traj.len\t"; // 4
    output << "cent.x\t"; // 5 
    output << "cent.y\t"; // 6 
    output << "cent.z\t"; // 7
    output << "surface.area\t"; // 8
    output << "nuc.apical.dist\t"; // 9
    output << "neighbor.trajectory.id\t"; // 10
    output << "neighbor.cent.x\t"; // 11
    output << "neighbor.cent.y\t"; // 12
    output << "neighbor.cent.z\t"; // 13
    output << "neighbor.nuc.apical.dist\t"; // 14
    output << "neighbor.surface.area\t"; // 15
    output << "neighbor.bend\t"; // 16
    output << "neighbor.contact.sa\t"; // 17
    output << "total.contact.sa\t"; // 18
    output << "recip.bend\t"; // 19
    output << "recip.sa" << endl; // 20

    for(int trajectory_id = 0; trajectory_id < (int)traj_lengths.size(); trajectory_id++) {
      if(traj_lengths[trajectory_id] <= 0) continue;
      for(size_t t = 0; t < stacks.size(); t++) {
	map<int,int> &traj_id2cell = stacks[t]->traj_id2cell;
	if(traj_id2cell.find(trajectory_id) == traj_id2cell.end()) continue;

	int cell_idx = traj_id2cell[trajectory_id];
	image_cell *cell = stacks[t]->cells.at(cell_idx);
	vec3 centroid = conf.voxel_alpha() * cell->centroid;
	
	vector<int> &neighbors = cell->neighbors;
	float total_contact_sa = 0;
	for(size_t n = 0; n < neighbors.size(); n++) total_contact_sa += cell->neighbor_contact_sa[n];

	for(size_t n = 0; n < neighbors.size(); n++) {
	  int neighbor_idx = neighbors[n];
	  image_cell *neighbor = stacks[t]->cells.at(neighbor_idx);
	  
	  output << t + conf.startframe + 1 << '\t'; // frame 1
	  output << (t + conf.startframe) * conf.time_step << '\t'; // time 2
	  output << trajectory_id << '\t'; // trajectory.id 3	  
	  output << traj_lengths[trajectory_id] << '\t'; // traj.len 4
	  output << centroid.x() << '\t'; // cent.x 5
	  output << centroid.y() << '\t'; // cent.y 6
	  output << centroid.z() << '\t'; // cent.z 7
	  output << cell->surface_area << '\t'; // surface.area 8
	  output << cell->nuc_apical_dist << '\t'; // nuc.apical.dist 9
	  output << neighbor->trajectoryid << '\t'; // neighbor.trajectory.id 10
	  vec3 neighbor_centroid = conf.voxel_alpha() * neighbor->centroid;
	  output << neighbor_centroid.x() << '\t'; // cent.x 11
	  output << neighbor_centroid.y() << '\t'; // cent.y 12
	  output << neighbor_centroid.z() << '\t'; // cent.z 13
	  output << neighbor->nuc_apical_dist << '\t'; // neighbor.nuc.apical.dist 14
	  output << neighbor->surface_area << '\t'; // neighbor.surface.area 15
	  output << cell->neighbor_bend[n] << '\t'; // neighbor.bend 16
	  output << cell->neighbor_contact_sa[n] << '\t';  // neighbor.contact.sa 17
	  output << total_contact_sa << '\t'; // 18

	  float recip_bend = -1, recip_sa = -1;
	  for(size_t r = 0; r < neighbor->neighbors.size(); r++) {
	    if(neighbor->neighbors[r] == cell_idx) {
	      if(neighbor->neighbor_bend.size() == neighbor->neighbors.size() &&
		 neighbor->neighbor_contact_sa.size() == neighbor->neighbors.size()) {
		recip_bend = neighbor->neighbor_bend[r];
		recip_sa = neighbor->neighbor_contact_sa[r];
	      }
	    }
	  }
	  output << recip_bend << '\t'; // recip.bend 19
	  output << recip_sa << endl; // recip.sa 20

	  // TODO: Scan neighbor->neighbors for cell_idx to get "recipricol contact"
	}
      }
    }
  }

  void hyperstack::save_analyze_trajectories(string filename) {
    measure_trajectories_threaded();

    ofstream output(filename.c_str());
    output << "frame\t"; // 1
    output << "time\t"; // 2 
    output << "trajectory.id\t";  // 3
    output << "traj.len\t"; // 4
    output << "cent.x\t"; // 5 
    output << "cent.y\t"; // 6 
    output << "cent.z\t"; // 7
    output << "volume\t"; // 8 
    output << "surface.area\t"; // 9 
    output << "AB.length\t"; // 10
    output << "anisotropy\t"; // 11A
    output << "anisotropy.AB\t"; // 11B    
    output << "vol.above.nuc\t"; // 12
    output << "sa.above.nuc\t"; // 13
    output << "apical.x\t"; // 14
    output << "apical.y\t"; // 15
    output << "apical.z\t"; // 16
    output << "basal.x\t"; // 17
    output << "basal.y\t"; // 18
    output << "basal.z\t"; // 19
    output << "cell.idx\t"; // 20
    output << "max.cell.bend\t"; // 21
    output << "nuc.apical.dist" << endl; // 22

    float alpha = conf.voxel_alpha();
    for(int trajectory_id = 0; trajectory_id < (int)traj_lengths.size(); trajectory_id++) {
      if(traj_lengths[trajectory_id] <= 0) continue;

      for(size_t t = 0; t < stacks.size(); t++) {
	map<int,int> &traj_id2cell = stacks[t]->traj_id2cell;
	if(traj_id2cell.find(trajectory_id) == traj_id2cell.end()) continue;

	int cell_idx = traj_id2cell[trajectory_id];
	image_cell *cell = stacks[t]->cells.at(cell_idx);
      
	vec3 centroid = alpha * cell->centroid;

	output << t + conf.startframe + 1 << '\t'; // frame 1
	output << (t + conf.startframe) * conf.time_step << '\t'; // time 2
	output << trajectory_id << '\t'; // trajectory.id 3
	output << traj_lengths[trajectory_id] << '\t'; // traj.len 4
	output << centroid.x() << '\t'; // cent.x 5
	output << centroid.y() << '\t'; // cent.y 6
	output << centroid.z() << '\t'; // cent.z 7
	output << cell->model_volume << '\t'; // volume 8
	output << cell->surface_area << '\t'; // surface.area
	vector<float> &pca_dims = cell->pca_dims;
	output << pca_dims.at(2) << '\t'; // AB.length 9
	float anisotropy = pca_dims.at(1) / pca_dims.at(0);
	output << anisotropy << '\t'; 	// anisotropy 11A
	float anisotropy_AB = pca_dims.at(2) / pca_dims.at(1);
	output << anisotropy_AB << '\t'; // anisotropy 11B
	output << cell->vol_above_nuc << '\t';	 // vol.above.nuc 12
	output << cell->sa_above_nuc << '\t'; // sa.above.nuc 13
	vec3 apical_pos = alpha * cell->apical_pos, basal_pos = alpha * cell->basal_pos;
	output << apical_pos.x() << '\t'; // 14
	output << apical_pos.y() << '\t'; // 15
	output << apical_pos.z() << '\t'; // 16
	output << basal_pos.x() << '\t'; // 17
	output << basal_pos.y() << '\t'; // 18
	output << basal_pos.z() << '\t'; // 19
	output << cell_idx << '\t'; // 20
	output << cell->max_cell_bend << '\t'; // 21
	output << cell->nuc_apical_dist << endl; // 22
      }        
    }
  }

  void image_stack::split_nucs() {
    cout << "running split_nucs()" << endl;
    vector< vector<int> > cells_assigned2nuc(nuc_voxels.size());
    
    volume32 nucs(iwidth, iheight, idepth);  nucs.seed(nuc_voxels);

    // Scan through voxels of each cell.
    for(int cell_idx = 0; cell_idx < (int)cell_voxels.size(); cell_idx++) {
      vector<ivec3> &cell_vox = cell_voxels[cell_idx];

      // Get overlapping nucs and keep track of the size of the overlap.
      map<int, int> nuc_idxs;
      for(size_t i = 0; i < cell_vox.size(); i++)  { 
	int nuc_idx = nucs(cell_vox[i]) - 1;
	if(nuc_idx >= 0) {
	  if( nuc_idxs.find(nuc_idx) != nuc_idxs.end() ) 
	    nuc_idxs[nuc_idx]++; 
	  else 
	    nuc_idxs[nuc_idx] = 1;
	}
      }

      // Count nucs that overlap > 10% of the nuc's volume.
      int overlapping_nuc_cnt = 0, onuc_idx = -1;
      for(map<int,int>::iterator it = nuc_idxs.begin(); it != nuc_idxs.end(); it++) {
	int nuc_idx = it->first;
	int nuc_overlap = it->second;
	vector<ivec3> &nuc_vox = nuc_voxels[nuc_idx];
	if(nuc_overlap > (int)nuc_vox.size() / 10) { overlapping_nuc_cnt++; onuc_idx = nuc_idx; }
      }

      // If there is only one such nuc, then, assign the cell to the nuc.
      if(overlapping_nuc_cnt == 1) cells_assigned2nuc[onuc_idx].push_back(cell_idx);
    }

    
    volume32 cells(iwidth, iheight, idepth);
    cells.seed(cell_voxels);

    // Set nuc regions to 0 and background to 255.
    volume8 fix_regions(iwidth, iheight, idepth);
    fix_regions.fill(255);
    for(size_t n = 0; n < nuc_voxels.size(); n++) {
      for(size_t i = 0; i < nuc_voxels[n].size(); i++) { fix_regions(nuc_voxels[n][i]) = 0; }
    }

    // Compute negated EDT in nuc regions only.
    volume32 fix_negEDT(iwidth, iheight, idepth);
    fix_negEDT.ComputeEDT(fix_regions);
    fix_negEDT.negate();
    
    // Mask out non-nuc regions in label set.
    volume32 fix_labels(iwidth, iheight, idepth);  
    fix_labels.seed(fix_regions, 255, -1);

    int fix_label = 1;
    for(int nuc_idx = 0; nuc_idx < (int)cells_assigned2nuc.size(); nuc_idx++) {
      vector<ivec3> &nuc_vox = nuc_voxels[nuc_idx];
      if(cells_assigned2nuc[nuc_idx].size() > 1) {
	// Split nucs that span multiple cells.
	for(size_t c = 0; c < cells_assigned2nuc[nuc_idx].size(); c++) {
	  // +1 here since labels are 1-indexed.
	  int cell_label = cells_assigned2nuc[nuc_idx][c] + 1;
	  int fix_label_cnt = 0;
	  // Look for regions in merged nuc that overlap this cell label.
	  for(size_t i = 0; i < nuc_vox.size(); i++) {
	    // Assign the new label.
	    if( cells(nuc_vox[i]) == cell_label ) { 
	      fix_labels(nuc_vox[i]) = fix_label;
	      fix_label_cnt++;
	    }
	  }
	  // Advance to next label.
	  fix_label++;
	}
      }
      else {
	// No assigned cells or a single cell assigned to this nuc, just apply label.
	for(size_t i = 0; i < nuc_vox.size(); i++) fix_labels(nuc_vox[i]) = fix_label;
	// Advance to next label.
	fix_label++;
      }
    }

    // Use Watershed to expand over small regions.
    fix_labels.watershed(fix_negEDT);

    // Get new set of nuclei.
    vector< vector<ivec3> > nuc_labels;
    fix_labels.components(nuc_labels);  // Find components with a label > 0.
    cout << "nuc_voxels.size() = " << nuc_voxels.size() << endl;
    nuc_voxels.clear(); 

    // Max nuc volume is only applied here.
    for(size_t n = 0; n < nuc_labels.size(); n++) { 
      if((int)nuc_labels[n].size() >= conf.min_nuc_vol && 
	 (int)nuc_labels[n].size() <= conf.max_nuc_vol) nuc_voxels.push_back(nuc_labels[n]);
    }
    cout << "nuc_labels.size() = " << nuc_voxels.size() << endl;

  }

  void image_stack::find_cell_neighbors2() {
    volume32 cell_labels(iwidth, iheight, idepth);
    // Note that these labels are 1-indexed, 0 = background
    cell_labels.seed(cell_voxels); 

    // Compute cell neighbors.
    cell_neighbors.resize(cell_voxels.size());
    cell_surface_voxels.resize(cell_voxels.size());
    cell_surface_neighbors.resize(cell_voxels.size());

    for(int cidx = 0; cidx < (int)cell_voxels.size(); cidx++) {
      cell_neighbors[cidx].clear(); // Clear any previous data.
      cell_surface_voxels[cidx].clear();
      cell_surface_neighbors[cidx].clear();

      // Get cell voxels.
      vector<ivec3> &cell_vox = cell_voxels[cidx];
      int cur_label = cidx + 1; // 1-indexed labels

      // Map from neighborighbor label to contact voxel count.
      map<int, int> nlabel2contact; 

      // Scan through voxels on cell
      for(size_t i = 0; i < cell_vox.size(); i++) {
	// Determine if the voxel is a surface voxel.
	vector<ivec3> hood; cell_labels.nhood26(hood, cell_vox[i]);
	
	bool is_surf_voxel = false;
	for(size_t h = 0; h < hood.size(); h++) {
	  if(cell_labels(hood[h]) == 0 || cell_labels(hood[h]) != cur_label) {
	    is_surf_voxel = true; break;
	  }
	}
	// Find which voxels should be checked.
	vector<ivec3> check;
	for(size_t h = 0; h < hood.size(); h++) {
	  // Check voxel since not current label and not zero.
	  if(cell_labels(hood[h]) != cur_label && cell_labels(hood[h]) != 0) check.push_back(hood[h]);
	  else if(cell_labels(hood[h]) == 0) {
	    // Thinned zero voxel. So expand region to check.
	    vector<ivec3> expand; 
	    cell_labels.nhood26(expand, hood[h]);
	    // Keep only voxels that are not zero and not current label.
	    for(size_t e = 0; e < expand.size(); e++) {
	      if(cell_labels(expand[e]) != cur_label && cell_labels(expand[e]) != 0) {
		bool found = false;
		for(size_t k = 0; k < check.size(); k++) {
		  if(check[k] == expand[e]) { found = true;  break; }
		}
		if(!found) check.push_back(expand[e]);
	      }
	    }
	  } 
	}
	if(check.size() == 0) {
	  if(is_surf_voxel) {
	    cell_surface_voxels[cidx].push_back(cell_vox[i]);
	    cell_surface_neighbors[cidx].push_back(-1);
	  }
	  continue; // Internal voxel, skip.
	}

	// Vote on a neighbor assignment.
	map<int, int> n2vote;
	for(size_t c = 0; c < check.size(); c++) {
	  int check_label = cell_labels(check[c]);
	  if(n2vote.find(check_label) == n2vote.end()) n2vote[check_label] = 1;
	  else n2vote[check_label]++;
	}

	int max_vote = 0, max_neighbor = 0;
	for(map<int,int>::iterator it = n2vote.begin(); it != n2vote.end(); it++) {
	  int neighbor = it->first;
	  int vote = it->second;
	  if(vote > max_vote) {
	    max_neighbor = neighbor;
	    max_vote = vote;
	  }
	}
	if(max_neighbor > 0) {
	  if(nlabel2contact.find(max_neighbor) == nlabel2contact.end()) 
	    nlabel2contact[max_neighbor] = 1;
	  else nlabel2contact[max_neighbor]++;
	}	
	cell_surface_voxels[cidx].push_back(cell_vox[i]);
	cell_surface_neighbors[cidx].push_back(max_neighbor-1);
      }      
      
      for(map<int,int>::iterator it = nlabel2contact.begin(); it != nlabel2contact.end(); it++) {
	int nlabel = it->first;
	int contact_count = it->second;
	// If contact > 1% of surface voxels, save as a nighbor.
	if(float(contact_count) > (float)cell_surface_voxels.size() * conf.neighbor_alpha) 
	  cell_neighbors[cidx].push_back(nlabel - 1);
      }

      /*
	volume8 surface(cell_surface_voxels[cidx]);
      ivec3 surf_shift(surface.x0, surface.y0, surface.z0);
      volume32 pos2idx(surface.width, surface.height, surface.depth);
      pos2idx.fill(-1);
      for(size_t n = 0; n < cell_surface_neighbors[cidx].size(); n++) {
	pos2idx(cell_surface_voxels[cidx][n] - surf_shift) = n;
      }
      vector<int> new_surf_neighbors = cell_surface_neighbors[cidx];
      for(size_t n = 0; n < cell_surface_neighbors[cidx].size(); n++) {
	if(cell_surface_neighbors[cidx][n] == -1) {
	  vector<ivec3> hood;
	  surface.nhood6(hood, cell_surface_voxels[cidx][n] - surf_shift);
	  for(size_t h = 0; h < hood.size(); h++) {
	    if(pos2idx(hood[h]) != -1) {
	      int idx = pos2idx(hood[h]);
	      new_surf_neighbors.at(idx) = -1;
	    }
	  }
	}
      }      
      cell_surface_neighbors[cidx] = new_surf_neighbors;*/

      // Clear neighbors that have spurious contact.
      for(size_t n = 0; n < cell_surface_neighbors[cidx].size(); n++) {
	bool found_neighbor = false;
	for(size_t i = 0; i < cell_neighbors[cidx].size(); i++) {
	  if(cell_surface_neighbors[cidx][n] == cell_neighbors[cidx][i]) { found_neighbor = true; break; }
	}
	if(found_neighbor == false) cell_surface_neighbors[cidx][n] = -1; // Clear this neighbor.
      }
      
    }

  }



  struct trajcontact_t {
    trajcontact_t() { trajid1 = trajid2 = timestep = -1;  }
    trajcontact_t(int tid1, int tid2, int time) {
      trajid1 = min(tid1, tid2);
      trajid2 = max(tid1, tid2);
      timestep = time;
    }
    int trajid1, trajid2, timestep;
    bool operator < (const trajcontact_t &B) const {
      if(trajid1 == B.trajid1) {
	if(trajid2 == B.trajid2) return timestep < B.timestep;
	else return trajid2 < B.trajid2;
      }
      else return trajid1 < B.trajid1;
    }
    bool operator == (const trajcontact_t &B) const {
      return trajid1 == B.trajid1 && trajid2 == B.trajid2 && timestep == B.timestep;
    }
  };

  void hyperstack::save_stable_neighbors(string filename) {
    ofstream output(filename.c_str());
    set<trajcontact_t> contacts;

    for(int trajectory_id = 0; trajectory_id < (int)traj_lengths.size(); trajectory_id++) {
      // Skip non-trajectories.
      if(traj_lengths[trajectory_id] <= 1) continue; 
      
      for(size_t t = 0; t < stacks.size(); t++) {
	// Find mapping from trajectory ID to cell.
	map<int,int> &traj_id2cell = stacks[t]->traj_id2cell;

	// Look for trajectory ID, skip timestep if none found.
	if(traj_id2cell.find(trajectory_id) == traj_id2cell.end()) continue;

	// Get cell.
	image_cell *cur_cell = stacks[t]->cells.at(traj_id2cell[trajectory_id]);

	// Save trajectory contacts of cell neighbors.
	vector<int> &neighbors = cur_cell->neighbors;
	for(size_t n = 0; n < neighbors.size(); n++) {
	  image_cell *neighbor = stacks[t]->cells.at(neighbors[n]);
	  if(neighbor->trajectoryid >= 0) {
	    // TODO: This needs to keep traj1 < traj 2.
	    trajcontact_t contact(trajectory_id, neighbor->trajectoryid, t);
	    contacts.insert(contact);
	  }
	}
      }
    }
    /*vector<trajcontact_t> contacts_vec(contacts.size());
      size_t vec_idx = 0;*/
    output << "traj.id1" << '\t' << "traj.id2" << '\t' << "len" << endl;
    int contact_len = 1;
    trajcontact_t prev_contact;
    bool first_flag = false;
    for(set<trajcontact_t>::iterator it = contacts.begin(); it != contacts.end(); it++) {
      trajcontact_t cur_contact = *it;
      if(!first_flag) { 
	first_flag = true; 
	prev_contact = cur_contact;
	continue; 
      }
      if(prev_contact.trajid1 == cur_contact.trajid1 && prev_contact.trajid2 == cur_contact.trajid2) {
	contact_len++;
      }
      else {
	output << prev_contact.trajid1 << '\t' << prev_contact.trajid2 << '\t' << contact_len << endl;
	contact_len = 1;
      }
      prev_contact = cur_contact;
    }
    output << prev_contact.trajid1 << '\t' << prev_contact.trajid2 << '\t' << contact_len << endl;
    
  }

  void hyperstack::save_neighbor_swaps(string filename) {
    ofstream output(filename.c_str());

    if(stacks.size() == 0) return;

    vector<int> swap_events(stacks.size(), 0), cells_analyzed(stacks.size(), 0);

    for(int trajectory_id = 0; trajectory_id < (int)traj_lengths.size(); trajectory_id++) {
      // Skip non-trajectories.
      if(traj_lengths[trajectory_id] <= 10) continue; 

      image_cell *cur_cell = NULL, *prev_cell = NULL;
      int prev_t = -1;
      for(size_t t = 0; t < stacks.size(); t++) {
	map<int,int> &traj_id2cell = stacks[t]->traj_id2cell;
	if(traj_id2cell.find(trajectory_id) == traj_id2cell.end())  continue;

	cur_cell = stacks[t]->cells.at(traj_id2cell[trajectory_id]);

	if(stacks[t]->is_border_touch(cur_cell->voxels)) break;

	if(prev_cell == NULL) { prev_cell = cur_cell; prev_t = t; continue;  }

	// Get trajectory IDs of neighbors for previous and current cell.
	vector<int> prev_trajids, cur_trajids;

	vector<int> &cur_neighbors = cur_cell->neighbors;
	for(size_t n = 0; n < cur_neighbors.size(); n++) {
	  image_cell *neighbor = stacks[t]->cells.at(cur_neighbors[n]);
	  if(neighbor->trajectoryid >= 0) cur_trajids.push_back(neighbor->trajectoryid);
	}

	vector<int> &prev_neighbors = prev_cell->neighbors;
	for(size_t n = 0; n < prev_neighbors.size(); n++) {
	  // Not quite right.
	  image_cell *neighbor = stacks[prev_t]->cells.at(prev_neighbors[n]);
	  if(neighbor->trajectoryid >= 0) prev_trajids.push_back(neighbor->trajectoryid);
	}

	int same_cur = 0, new_neighbors = 0;
	for(size_t i = 0; i < cur_trajids.size(); i++) {
	  bool found_match = false;
	  for(size_t j = 0; j < prev_trajids.size(); j++) {
	    if(cur_trajids[i] == prev_trajids[j]) {
	      found_match = true;
	      break;
	    }
	  }
	  if(found_match) same_cur++; else new_neighbors++;
	}

	int same_prev = 0, lost_neighbors = 0;
	for(size_t i = 0; i < prev_trajids.size(); i++) {
	  bool found_match = false;
	  for(size_t j = 0; j < cur_trajids.size(); j++) {
	    if(prev_trajids[i] == cur_trajids[j]) {
	      found_match = true;
	      break;
	    }
	  }
	  if(found_match) same_prev++; else lost_neighbors++;
	}
	// More than half the neighbors changed, count as a swap event.
	if(float(lost_neighbors) / float(same_cur + lost_neighbors) > 0.5) {
	  swap_events[t]++;
	}
	cells_analyzed[t]++;

	prev_cell = cur_cell;
	prev_t = t;
      }
    }    
    
    output << "frame\t" << "time\t" << "cells.analyzed\t" << "swap.events" << endl;
    for(size_t t = 1; t < stacks.size(); t++) {
      output << t + conf.startframe + 1 << '\t';
      output << (t + conf.startframe) * conf.time_step << '\t';
      output << cells_analyzed[t] << '\t';
      output << swap_events[t] << endl;
    }

  }

  void hyperstack::save_neighbor_analysis(string filename) {
    ofstream output(filename.c_str());

    output << "trajectory.id\t";
    for(size_t i = 0; i < 12; i++) output << "rank.contact." << i + 1 << '\t';
    output << "traj.len" << endl;

    // Scan through trajectories/cells.
    for(int trajectory_id = 0; trajectory_id < (int)traj_lengths.size(); trajectory_id++) {
      // Skip short trajectories.
      if(traj_lengths[trajectory_id] <= 1) continue; 

      int traj_len = 0;

      // Determine number of contacts with a trajectory.
      map<int, int> n_traj2cnt;
      for(size_t t = 0; t < stacks.size(); t++) {
	// Find mapping from trajectory ID to cell.
	map<int,int> &traj_id2cell = stacks[t]->traj_id2cell;

	// Look for trajectory ID, skip timestep if none found.
	if(traj_id2cell.find(trajectory_id) == traj_id2cell.end()) continue;

	// Get cell.
	image_cell *cur_cell = stacks[t]->cells.at(traj_id2cell[trajectory_id]);
	
	// Keep border contacts. if(stacks[t]->is_border_touch(cur_cell->voxels)) continue;
	traj_len++;

	// Save trajectory contacts of cell neighbors.
	vector<int> &neighbors = cur_cell->neighbors;
	for(size_t n = 0; n < neighbors.size(); n++) {
	  image_cell *neighbor = stacks[t]->cells.at(neighbors[n]);
	  if(neighbor->trajectoryid >= 0) {
	    // Update contact duration.
	    if(n_traj2cnt.find(neighbor->trajectoryid) == n_traj2cnt.end()) n_traj2cnt[neighbor->trajectoryid] = 1;
	    else n_traj2cnt[neighbor->trajectoryid]++;
	  }
	}
      }

      // Save and sort contacts along trajectory.
      vector<int> contact_cnts;
      for(map<int, int>::iterator it = n_traj2cnt.begin(); it != n_traj2cnt.end(); it++) 
	contact_cnts.push_back(it->second);

      sort(contact_cnts.begin(), contact_cnts.end(), greater<int>());

      output << trajectory_id << '\t';
      for(size_t i = 0; i < 12; i++) {
	if(contact_cnts.size() >= i+1) output << contact_cnts[i] << '\t';  else output << "NA\t";
      }
      //output << traj_lengths[trajectory_id] << endl;
      output << traj_len << endl;
    }
    // ------------
  }


  void hyperstack::save_neighbor_counts(string filename) {
    ofstream output(filename.c_str());

    output << "frame\t" << "time\t" << "cell.id\t" << "neighbor.cnt" << endl;

    for(size_t t = 0; t < stacks.size(); t++) {
      vector<image_cell *> &cells = stacks[t]->cells;

      for(size_t c = 0; c < cells.size(); c++) {
	if(cells[c]->trajectoryid < 0) continue;
	if(stacks[t]->is_border_touch(cells[c]->voxels)) continue;

	output << t + conf.startframe + 1 << '\t'; // frame
	output << (t + conf.startframe) * conf.time_step << '\t'; // time
	output << c << '\t'; // cell.id
	output << cells[c]->neighbors.size() << endl;
      }
    }

  }

  void image_stack::get_neighbors_traj(vector<int> &neighbor_traj, neighbormode_t neighbor_mode, 
				       vector<int> &orders, vector<int> &cur_traj) {
    if(cur_traj.size() == 0) return;
    vector<int> cur_idx, neighbor_idx;

    // Map from current trajectory to cell index.
    for(size_t c = 0; c < cur_traj.size(); c++) {
      if(traj_id2cell.find(cur_traj[c]) != traj_id2cell.end()) cur_idx.push_back(traj_id2cell[cur_traj[c]]);
    }    
    // Get order-neighbors using cell indices.
    get_neighbors(neighbor_idx, neighbor_mode, orders, cur_idx);

    if(traj_ids.size() > 0 && traj_ids.size() == cell_voxels.size()) {
      for(size_t n = 0; n < neighbor_idx.size(); n++) {
	int neighbor_traj_id = traj_ids[neighbor_idx[n]];
	if(neighbor_traj_id >= 0)
	  neighbor_traj.push_back(neighbor_traj_id);      
      }
    }
  }

  // Performs BFS to get n-order neighbors. 
  void image_stack::get_neighbors(vector<int> &neighbor_idx, neighbormode_t neighbor_mode, vector<int> &orders, vector<int> &cur_idx) {
    if(cell_voxels.size() == 0) return; 
    if(cur_idx.size() == 0) return;

    vector<vec3> cell_centroids(cell_voxels.size()); 
    for(size_t i = 0; i < cell_voxels.size(); i++) cell_centroids[i] = geom::Centroid3(cell_voxels[i]);

    vector<int> bfs_labels(cell_voxels.size(), -1);
    queue<int> Q;

    for(size_t c = 0; c < cur_idx.size(); c++) {
      bfs_labels.at(cur_idx[c]) = 0;
      Q.push(cur_idx[c]);
    }

    while(!Q.empty()) {
      int cur_idx = Q.front(); 
      Q.pop();
      vector<int> &neighbors = cell_neighbors[cur_idx];
      int cur_label = bfs_labels[cur_idx];
      for(size_t n = 0; n < neighbors.size(); n++) {
	if(bfs_labels[neighbors[n]] >= 0) continue;

	bool is_neighbor = false;
	switch(neighbor_mode) {
	case Neighbor_All: is_neighbor = true; break;
	case Neighbor_Xgreater: is_neighbor = cell_centroids[cur_idx].x() < cell_centroids[neighbors[n]].x(); break;
	case Neighbor_Xless: is_neighbor = cell_centroids[cur_idx].x() > cell_centroids[neighbors[n]].x(); break;
	default: cerr << "Unspported neighbor_mode" << endl; exit(1); 
	}
	if(is_neighbor ) {
	  bfs_labels[neighbors[n]] = cur_label + 1;
	  Q.push(neighbors[n]);
	}
      }
      for(size_t c = 0; c < cell_voxels.size(); c++) cells.at(c)->bfs_label = bfs_labels[c];
    }    
    
    for(int cell_idx = 0; cell_idx < (int)bfs_labels.size(); cell_idx++) {
      bool found_order = false;
      for(size_t o = 0; o < orders.size(); o++) { if(bfs_labels[cell_idx] == orders[o]) found_order = true; }
      if(found_order) 
	neighbor_idx.push_back(cell_idx);
    }
  }

  void hyperstack::save_nuc_cents(const string &filename) {
    ofstream output(filename.c_str());
    output << "frame" << '\t';
    output << "x" << '\t';
    output << "y" << '\t';
    output << "z" << endl;
    for(size_t t = 0; t < stacks.size(); t++) {
      int frame = stacks[t]->frame_num;
      vector<vec3> cents;
      stacks[t]->nuc_centroids(cents);
      for(size_t c = 0; c < cents.size(); c++) {
	output << frame << '\t';
	output << cents[c][0] << '\t' << cents[c][1] << '\t' << cents[c][2] << endl;
      }
    }    
  }
  
  void hyperstack::save_dorsalfolds(string filename, vector<int> &traj_ids) {
    ofstream output(filename.c_str());
    output << "frame" << '\t'; // -1
    output << "trajectory.id" << '\t'; // 0
    output << "id" << '\t'; // 1
    output << "total.surface.area" << '\t'; // 2
    output << "total.volume" << '\t'; // 3
    output << "above.baz.surface.area" << '\t'; //4 
    output << "above.baz.volume" << '\t'; // 5
    output << "cent.x" << '\t'; //6 
    output << "cent.y" << '\t'; //7
    output << "cent.z" << '\t'; //8
    output << "total.length" << '\t'; //9 
    output << "apical.length" << '\t'; //10
    output << "bazooka.pos" << '\t'; // 11
    output << "total.par.intensity" << '\t'; // 12
    output << "basal.par.intensity" << '\t'; // 13
    output << "bazooka.patch.volume" << '\t'; // 14
    output << "bazooka.patch.intensity" << '\t'; // 15
    output << "bazooka.nonpatch.volume" << '\t'; // 16
    output << "bazooka.nonpatch.intensity" << '\t'; // 17
    output << "basal.4voxel.par.volume" << '\t'; // 18
    output << "basal.4voxel.par.intensity" << '\t'; // 19
    output << "basal.2voxel.par.volume" << '\t'; // 20
    output << "basal.2voxel.par.intensity" << endl; // 21

    for(size_t t = 0; t < stacks.size(); t++) {
      map<int,int> &traj_id2cell = stacks[t]->traj_id2cell;

      for(size_t i = 0; i < traj_ids.size(); i++) {
	if(traj_id2cell.find(traj_ids[i]) == traj_id2cell.end()) continue;

	int cell_idx = traj_id2cell[traj_ids[i]];
	image_cell *cell = stacks[t]->cells.at(cell_idx);

	output << t + conf.startframe + 1 << '\t'; // frame -1
	output << cell->trajectoryid << '\t'; // 0
	output << cell->idx + 1 << '\t'; // 1
	output << cell->dorsalfolds.total_surface_area << '\t'; //2 
	output << cell->dorsalfolds.total_volume << '\t';  //3 
	output << cell->dorsalfolds.above_baz_surface_area << '\t'; // 4
	output << cell->dorsalfolds.above_baz_volume << '\t'; // 5
	output << cell->dorsalfolds.cent_x << '\t'; //6
	output << cell->dorsalfolds.cent_y << '\t'; //7
	output << cell->dorsalfolds.cent_z << '\t'; //8
	output << cell->dorsalfolds.total_length << '\t'; //9 
	output << cell->dorsalfolds.apical_length << '\t'; //10
	output << cell->dorsalfolds.bazooka_pos << '\t'; // 11
	output << fixed << cell->dorsalfolds.total_par_intensity << '\t'; // 12
	output << fixed << cell->dorsalfolds.basal_par_intensity << '\t'; // 13
	output << cell->dorsalfolds.bazooka_patch_volume << '\t'; // 14
	output << fixed << cell->dorsalfolds.bazooka_patch_intensity << '\t'; // 15
	output << cell->dorsalfolds.bazooka_nonpatch_volume << '\t'; // 16
	output << fixed << cell->dorsalfolds.bazooka_nonpatch_intensity << '\t'; // 17
	output << cell->dorsalfolds.basal_4voxel_par_volume << '\t'; // 18
	output << fixed << cell->dorsalfolds.basal_4voxel_par_intensity << '\t'; // 19
	output << cell->dorsalfolds.basal_2voxel_par_volume << '\t'; // 20
	output << fixed << cell->dorsalfolds.basal_2voxel_par_intensity << endl; // 21
      }
    }
  }

}


