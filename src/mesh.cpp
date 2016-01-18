#include <stdlib.h>
#include <limits>
#include <map>
#include <queue>
#include <set>

#include "mesh.hpp"
#include "util.hpp"

namespace geom {
  // Intersect two face lists. Return face index that does not match given face index.
  inline int isect(int fidx, vector<int> &A, vector<int> &B) {
    for(size_t i = 0; i < A.size(); i++) 
      for(size_t j = 0; j < B.size(); j++) if(A[i] == B[j] && A[i] != fidx) return A[i]; 
    // Nothing found in intersection, return -1.
    return -1;
  }

  // Compute face-face adjacencies by intersecting face lists on each vertex.
  void mesh::face_adjacencies() {
    for(size_t f = 0; f < faces.size(); f++) {
      int v0 = faces[f].vert[0], v1 = faces[f].vert[1], v2 = faces[f].vert[2];
      faces[f].adj[0] = isect(f, facelist[v0], facelist[v1]);
      faces[f].adj[1] = isect(f, facelist[v0], facelist[v2]);
      faces[f].adj[2] = isect(f, facelist[v1], facelist[v2]);
    }
  }

  // Update the mesh bounding box.
  void mesh::reset_range() {
    x_min = numeric_limits<float>().max(); x_max = -(numeric_limits<float>().max());
    y_min = numeric_limits<float>().max(); y_max = -(numeric_limits<float>().max());
    z_min = numeric_limits<float>().max(); z_max = -(numeric_limits<float>().max());

    for(size_t vidx = 0; vidx < vertices.size(); vidx++) {
      vec3 &v = vertices[vidx].pos;
      // Collect bounding box for mesh.
      if(v[0] > x_max) x_max = v[0]; if(v[0] < x_min) x_min = v[0];
      if(v[1] > y_max) y_max = v[1]; if(v[1] < y_min) y_min = v[1];
      if(v[2] > z_max) z_max = v[2]; if(v[2] < z_min) z_min = v[2];
    }
    // Set the center of the mesh
    x_center = x_min + (x_max - x_min) / 2; y_center = y_min + (y_max - y_min) / 2;  z_center = z_min + (z_max - z_min) / 2;

    centroid = compute_centroid();
    // Populate verticies.
    vector<vec3> vec3vert(vertices.size());
    for(size_t vidx = 0; vidx < vertices.size(); vidx++) vec3vert[vidx] = vertices[vidx].pos;

    PCA3(pca3, pca3_x, vec3vert); // Compute PCA of vertices.
    compute_axes();
  }

  void mesh::compute_axes() {
    // Populate verticies.
    vector<vec3> vec3vert(vertices.size());
    for(size_t vidx = 0; vidx < vertices.size(); vidx++) vec3vert[vidx] = vertices[vidx].pos - centroid;
    // Find maximum projection onto PCA axes.
    axes3.resize(3); 
    for(size_t i = 0; i < pca3.size(); i++) axes3[i] = pca3[i].max_mag_project(vec3vert);
  }

  void mesh::rebuild_adj() {
    // Compute how much space is needed for the face list on each vertex.
    vector<int> fadjcnt(vertices.size(), 0);
    for(size_t fidx = 0; fidx < faces.size(); fidx++) {
      for(int i = 0; i < 3; i++) fadjcnt[faces[fidx].vert[i]]++;
    }
    for(size_t vidx = 0; vidx < vertices.size(); vidx++) 
      if(fadjcnt[vidx] == 0) { cerr << "no adjacent face for vertex" << endl; exit(1); }

    // Reserve space for each face list.
    facelist.resize(vertices.size());

    // Clear, if any are present.
    for(size_t vidx = 0; vidx < facelist.size(); vidx++) facelist[vidx].clear();
    for(size_t vidx = 0; vidx < facelist.size(); vidx++) facelist[vidx].reserve(fadjcnt[vidx]); 

     // Populate each face list.
     for(int fidx = 0; fidx < (int)faces.size(); fidx++) {
       for(int i = 0; i < 3; i++) facelist[faces[fidx].vert[i]].push_back(fidx);
     }

    // Intersect face lists to obtain face-face adjacencies.
    face_adjacencies();  
  }

  // Create a mesh from a vertex table and a face table.
  void mesh::rebuild(vector<vec3> &vtable, vector<face> &ftable) {
    faces.clear(); vertices.clear();
    
    // Copy vertex positions.
    vertices.resize(vtable.size());
    // Copy mesh vertices.
    for(size_t vidx = 0; vidx < vtable.size(); vidx++) vertices[vidx].pos = vtable[vidx];

    // Copy verticies from face table.
    faces.resize(ftable.size());
    for(int fidx = 0; fidx < (int)faces.size(); fidx++) {
      for(int i = 0; i < 3; i++) faces[fidx].vert[i] = ftable[fidx].vert[i];
    }
    
    rebuild_adj(); // (re)Build face adjacencies
    update_meta(); // Update meta information.
  }

  // Compute face normals using vertex indicies.
  void mesh::compute_face_normals() {
    for(int fidx = 0; fidx < (int)faces.size(); fidx++) {
      mesh_face &F = faces[fidx];
      vec3 &v0 = vertices[F.vert[0]].pos;
      vec3 &v1 = vertices[F.vert[1]].pos;
      vec3 &v2 = vertices[F.vert[2]].pos;
      faces[fidx].normal = ComputeNormal(v0, v1, v2);
    }
  }

  // Computes vertex normal as the mean of the neighboring face normals.
  void mesh::compute_vertex_normals() {
    // Compute vertex normals. 
    for(size_t vidx = 0; vidx < vertices.size(); vidx++) {
      vector<int> &adj = facelist[vidx]; // Get face adjacencies.
      vec3 vnormal(0,0,0); 
      for(size_t a = 0; a < adj.size(); a++) vnormal += faces[adj[a]].normal;
      float Z = (float)adj.size();
      vertices[vidx].normal = (1.0f / Z) * vnormal;
      vertices[vidx].normal.normalize();
    }
  }

  // Voxelize mesh.
  void mesh::volpts(vector<ivec3> &pts) {
    // Create a binary volume using mesh voxelization.
    int w = ceil(x_max - x_min), h = ceil(y_max - y_min), d = ceil(z_max - z_min);
    ivec3 min_v(x_min, y_min, z_min);
    volume8 vol(w,h,d); voxelize(vol, 0, 255, true);

    for(int z = 0; z < d; z++)
      for(int y = 0; y < h; y++)
	for(int x = 0; x < w; x++) 
	  if(vol(x,y,z) == 255) pts.push_back(ivec3(x,y,z) + min_v);
  }

  // Voxalize, but only a thin boundary along edge of mesh.
  void mesh::volpts(vector<ivec3> &pts, int dthresh) {
    // Create a binary volume using mesh voxelization.
    int w = ceil(x_max - x_min), h = ceil(y_max - y_min), d = ceil(z_max - z_min);
    ivec3 min_v(x_min, y_min, z_min);
    // Voxelize
    volume8 vol(w,h,d); voxelize(vol, 0, 255, true);

    // Compute EDT, which allows selecting points along thin boundary.
    volume32 EDT(w,h,d); EDT.ComputeEDT(vol, 0);
    int dsqthresh = dthresh * dthresh;
    for(int z = 0; z < d; z++)
      for(int y = 0; y < h; y++)
	for(int x = 0; x < w; x++) {
	  if(vol(x,y,z) == 255) {
	    if(EDT(x,y,z) <= dsqthresh) {
	      pts.push_back(ivec3(x,y,z) + min_v);
	    }
	  }
	}
  }

  void mesh::volpts(vector<ivec3> &above, vector<ivec3> &below, plane3 &plane) {
    vector<ivec3> pts; volpts(pts);
    for(size_t p = 0; p < pts.size(); p++) {
      vec3 pt_v = pts[p].to_vec3();
      if(plane.pos_side(pt_v)) above.push_back(pts[p]);
      else if(plane.neg_side(pt_v)) below.push_back(pts[p]);
    }
  }

  double mesh::compute_volume(double vxd) {
    /// Make sure the mesh is closed.
    for(size_t fidx = 0; fidx < faces.size(); fidx++) {
      if(faces[fidx].adj[0] < 0) return -1; if(faces[fidx].adj[1] < 0) return -1;
      if(faces[fidx].adj[2] < 0) return -1;
    }
    // Create a binary volume for measurement.
    int w = ceil(x_max - x_min), h = ceil(y_max - y_min), d = ceil(z_max - z_min);
    // Voxelize the mesh in that volume.
    volume8 vol(w,h,d); voxelize(vol, 0, 255, true);
    double vxd3 = vxd * vxd * vxd, total_vol = 0;
    int vcnt = 0;
    // Add up the fg voxels in that small volume.
    for(int z = 0; z < d; z++)
      for(int y = 0; y < h; y++)
	for(int x = 0; x < w; x++)
	  if(vol(x,y,z) == 255) {
	    total_vol += vxd3;
	    vcnt++;
	  }
    
    //cout << "vol = " << vxd << "/" << vcnt << "/" << w * h * d << endl;
    return total_vol;
  }

  // Find areas of labeled components on the mesh.
  void mesh::component_areas(vector<int> &facelabels, float alpha, vector<float> &areas) {
    int num_components = 0;
    for(size_t i = 0; i < facelabels.size(); i++) num_components = max(num_components, facelabels[i]);
    
    areas.resize(num_components + 1);
    for(size_t a = 0; a < areas.size(); a++) areas[a] = 0;

    for(size_t fidx = 0; fidx < facelabels.size(); fidx++) {
      if(facelabels[fidx] >= 0) {
	mesh_face &F = faces[fidx];
	// This maps from voxel space to geometric space.
	vec3 p0 = alpha * vertices[F.vert[0]].pos;
	vec3 p1 = alpha * vertices[F.vert[1]].pos;
	vec3 p2 = alpha * vertices[F.vert[2]].pos;
	areas.at(facelabels[fidx]) += geom::area(p0, p1, p2);
      }
    }
  }

  // Clean components on surface image.
  void mesh::clean_small_comps(vector<int> &facelabels, int num_components, int thresh) {
    // Count the number of faces in each component.
    vector<int> labelcnt(num_components, 0);
    for(size_t i = 0; i < facelabels.size(); i++) if(facelabels[i] > 0) labelcnt.at(facelabels[i] - 1)++;

    // For those with less than threshold number of faces, set
    // facelabel to zero.
    vector< vector<int> > comps(num_components);
    for(size_t i = 0; i < facelabels.size(); i++) {
      if(facelabels[i] > 0 && labelcnt[facelabels[i] - 1] < thresh) { 
	facelabels[i] = 0;
      }
    }
  }

  // input labels, 0 = background and do not label 
  // 1 = label this triangle component
  // output, 0 = background and 1 - N component labels.
  int mesh::label_components(vector<int> &labels) {
    int num_components = 0;
    for(size_t fidx = 0;  fidx < labels.size(); fidx++) {
      if(labels[fidx] == 1) {
	dfs_label(labels, fidx, num_components+2);
	num_components++;
      }
    }
    // Subtract 1 so that labels range from [1,N].
    for(size_t i = 0; i < labels.size(); i++) if(labels[i] != 0) labels[i] -= 1;
    return num_components;
  }

  void mesh::connected_components(vector< vector<int> > &components, int min_tri) {
    // Create a label for each face.
    vector<int> labels(faces.size(), 1); // use 1 here to ensure all faces visited
    int num_components = label_components(labels);
    vector< vector<int> > components_raw(num_components);
    for(size_t i = 0; i < labels.size(); i++) {
      if(labels[i] > 0) { components_raw.at(labels[i] - 1).push_back(i); }
    }
    for(size_t i = 0; i < components_raw.size(); i++) {
      if((int)components_raw[i].size() > min_tri) components.push_back(components_raw[i]);
    }
  }

  // Label nodes in a connected component by DFS.  0 = background, 1 =
  // to label, SCANNED = in queue.
  void mesh::dfs_label(vector<int> &labels, int start, int label) {
    const int SCANNED = -numeric_limits<int>::max();
    vector<int> L; L.push_back(start); 
    labels[start] = SCANNED; // SCANNED = in queue
    while(L.empty() == false) {
      // Get current face from queue.
      int cur = L.back(); L.pop_back();
      labels[cur] = label; // mark as finished

      // Get adjacent faces.
      int *adj = faces[cur].adj;
      for(int j = 0; j < 3; j++) {
	if(adj[j] < 0) continue; // Skip edges that have no adjacencies.
	if(labels[adj[j]] == 1) {  // unvisited 
	  labels[adj[j]] = SCANNED; // SCANNED means it's in the queue
	  L.push_back(adj[j]); 
	}
      }
    }
  }

  // Starting a given face index start, label face with labe with a
  // given threshold distance.
  void mesh::dfs_paint(vector<int> &facelabels, int start, int label, int dthresh) {
    vector<int> L, dist; L.push_back(start); dist.push_back(0);
    // Uses most negative value for labeling faces that got scanned.
    const int SCANNED = -numeric_limits<int>::max(), VISITED = -numeric_limits<int>::max() + 1;
    facelabels[start] = SCANNED;

    vector<int> to_label;
    while(L.empty() == false) {
      // Get current face from queue.
      int cur = L.back(); L.pop_back();
      // triangle cur_tri = as_triangle(cur);
      // vec3 cur_cent = cur_tri.centroid();
      int cur_dist = dist.back(); dist.pop_back();
      facelabels[cur] = VISITED;
      to_label.push_back(cur);

      if(cur_dist < dthresh) {
	// Get adjacent faces.
	int *adj = faces[cur].adj;
	for(int j = 0; j < 3; j++) {
	  if(adj[j] < 0) continue; // Skip edges that have no adjacencies.
	  if(facelabels[adj[j]] == SCANNED) continue; // In queue, skip.
	  if(facelabels[adj[j]] == VISITED) continue; // Already visited, skip.
	  // triangle adj_tri = as_triangle(adj[j]);
	  // float d = geom::distance3(adj_tri.centroid(), cur_cent);
	  dist.push_back(cur_dist + 1);
	  L.push_back(adj[j]); 
	}
      }
    }
    // Apply label to visited faces.
    for(size_t i = 0; i < to_label.size(); i++) facelabels[to_label[i]] = label;
  }

  // Hash map that implements linear probing.  TODO: Make this a hash table that grows.
  struct entry { int vidx_mesh, vidx_sub; }; 
  struct hash_map {
    int N; entry *data;
    hash_map(int Nmax) { N = Nmax; data = new entry[N]; clear(); }
    ~hash_map() { delete[] data; }
    // if vidx_mesh = -1, then the slot is empty.
    void clear() { for(int i = 0; i < N; i++) data[i].vidx_mesh = -1; }
    int find(int vidx_mesh) {
      int i = vidx_mesh % N;
      int istart = i;
      while(data[i].vidx_mesh >= 0) { 
	if(data[i].vidx_mesh == vidx_mesh) return data[i].vidx_sub; 
	else i = (i + 1) % N;
	if(i == istart) return -1; 
      }
      return -1;
    }
    void insert(int vidx_mesh, int vidx_sub) {
      int i = vidx_mesh % N;
      int istart = i; 
      while(data[i].vidx_mesh >= 0) { 
	i = (i + 1) % N; 
	if(i == istart) { cerr << "mesh: hashmap_out of space" << endl; exit(1); }
      }
      data[i].vidx_mesh = vidx_mesh; data[i].vidx_sub = vidx_sub;
    }
  };

  // Returns a component as a mesh face table and coresponding vertex table.
  void mesh::as_mesh(vector<int> &subfaces, vector<vec3> &vtable, vector<face> &ftable) {
    vtable.clear(); ftable.clear();
    hash_map vidx(3*subfaces.size());
    for(size_t i = 0; i < subfaces.size(); i++) {
      // Get face index.
      int fidx = subfaces[i]; 
      // Get vertex indicies of that face.
      int v0 = faces[fidx].vert[0], v1 = faces[fidx].vert[1], v2 = faces[fidx].vert[2];

      // Map those vertex indicies to vertiex indicies in the sub-mesh.
      int v0sub = vidx.find(v0), v1sub = vidx.find(v1), v2sub = vidx.find(v2);
      // If any of the verticies were not found in the submesh, create
      // a new vertex in the submesh and save the map from the old
      // vertex to the new one.
      if(v0sub < 0) {
	v0sub = vtable.size();
	vtable.push_back(vertices[v0].pos);
	vidx.insert(v0, v0sub);
      }
      if(v1sub < 0) {
	v1sub = vtable.size();
	vtable.push_back(vertices[v1].pos);
	vidx.insert(v1, v1sub);
      }
      if(v2sub < 0) {
	v2sub = vtable.size();
	vtable.push_back(vertices[v2].pos);
	vidx.insert(v2, v2sub);
      }
      // Create new face with vertex indicies in the submesh and copy the normal.
      face newface;
      // newface.normal = faces[fidx].normal;
      newface.vert[0] = v0sub; newface.vert[1] = v1sub; newface.vert[2] = v2sub;
      ftable.push_back(newface);
    }
  }


  vec3 mesh::refine_disp(volume8 &v, int vidx, float dist, float sign) {
    vec3 o = vertices[vidx].pos;
    vec3 d = vertices[vidx].normal;
    d = sign * d; // Normals point inward, this is weird so flip sign.
    vec3 p = o; // Initially no displacement.
    // Intersect volume within some threshold distance.
    vector<ivec3> res; v.ray_intersect(res, o, d, dist); 
    if(res.size() > 0) {
      // If items found, then find maximum intensity voxel.
      uint8 maxI = v(res[0]); ivec3 maxpos = res[0];
      for(size_t r = 1; r < res.size(); r++) {
	if( v(res[r]) > maxI ) { maxI = v(res[r]); maxpos = res[r]; }
      }
      p[0] = float(maxpos.x()); p[1] = float(maxpos.y()); p[2] = float(maxpos.z());
    }
    p -= o;
    return p; // Compute displacement to position.
  }

  // TODO: Use smooth vertex and maxI vertex add them together to get
  // new position add alpha/(1-alpha) weight to smoothing vs fitting.
  // Refines mesh position to match high intensity image boundaries in
  // source volume v.
  void mesh::refine(volume8 &v, float dist, int max_iter, float alpha, float stepsize) {
    // Vertex displacements.
    vector<vec3> sdisp(vertices.size()); 
    vector<vec3> disp(vertices.size()); 
    for(int iter = 0; iter < max_iter; iter++) {
      for(int vidx = 0; vidx < num_vert(); vidx++) {
	disp[vidx] = refine_disp(v, vidx, dist);
	sdisp[vidx] = smooth_disp(vidx);
      }
      for(int vidx = 0; vidx < num_vert(); vidx++) 
	vertices[vidx].pos += stepsize * (alpha * disp[vidx] + (1.0 - alpha) * sdisp[vidx]);
      // Vertices moved so re-compute face normals.
      compute_face_normals(); compute_vertex_normals();
    }
    update_meta();
  }

  vec3 mesh::smooth_disp(int vidx) {
    vec3 &v = vertices[vidx].pos; // v = current vertex 
    vector<int> &fadj = facelist[vidx];  // get face adjacencies for the current vertex
    int Z = 0; 
    vec3 smoothed(0,0,0);
    // Get faces adjacent to that vertex.
    for(int f = 0; f < (int)fadj.size(); f++) {
      mesh_face &adjface = faces[fadj[f]];
      // Examine the vertices of that face.
      for(int i = 0; i < 3; i++) {
	if(adjface.vert[i] != vidx) { 
	  smoothed += vertices[adjface.vert[i]].pos;
	  Z++;
	}
      }
    }
    if(Z == 0) return v;
    smoothed /= (float)Z;
    smoothed -= v;
    return smoothed;
  }

  void mesh::smooth_iter(vector<vec3> &disp, float stepsize) {
    if(disp.size() != vertices.size()) disp.resize(vertices.size());
    for(int vidx = 0; vidx < num_vert(); vidx++) {
      vec3 &v = vertices[vidx].pos; // v = current vertex 
      vector<int> &fadj = facelist[vidx];  // get face adjacencies for the current vertex
      int Z = 0; 
      vec3 smoothed(0,0,0);
      // Get faces adjacent to that vertex.
      for(int f = 0; f < (int)fadj.size(); f++) {
	mesh_face &adjface = faces[fadj[f]];
	// Examine the vertices of that face.
	for(int i = 0; i < 3; i++) {
	  if(adjface.vert[i] != vidx) { 
	    smoothed += vertices[adjface.vert[i]].pos;
	    Z++;
	  }
	}
      }
      if(Z == 0) continue;
      // The smoothed position is just the average of those face verticies.
      smoothed /= (float)Z;
      smoothed -= v;
      disp[vidx] = smoothed; // Compute and store the displacement to smoothed position.
    }
    // Displace face in direction of smoothed position by given step size. 
    for(int vidx = 0; vidx < num_vert(); vidx++) vertices[vidx].pos += stepsize * disp[vidx];     
    // Recompute vertex and face normals and advance to next iteration.
    compute_face_normals();
  }

  // Performs Taubin lambda/mu mesh smoothing. Algorithm adapted from
  // TriMesh2.
  void mesh::smooth(int niters) {
    float stepsize = 0;
    vector<vec3> disp(vertices.size()); // Vertex displacements.
    for(int iter = 0; iter < 2 * niters; iter++) {
      if(iter % 2 == 0) stepsize = 0.330f; else stepsize = -0.331f;
      smooth_iter(disp, stepsize);
    }
    update_meta();
  }

  // For each face create a triangle that has a map from the triangle
  // to the face index. Build a BVH on those triangles.
  void mesh::build_bvh() {
    if(bvh != NULL) delete bvh;
    tri_face.resize(faces.size());
    for(int fidx = 0; fidx < (int)tri_face.size(); fidx++) {
      mesh_face &face = faces[fidx];
      mesh_triangle t;
      t.normal = face.normal;
      t.v[0] = vertices[face.vert[0]].pos; t.v[1] = vertices[face.vert[1]].pos; t.v[2] = vertices[face.vert[2]].pos;
      t.fidx = fidx;
      tri_face[fidx] = t;
    }
    bvh = new BVH<mesh_triangle>(tri_face);
  }

  // Performs ray intersection using the BVH. Returns intersection
  // points and face indexes.
  void mesh::ray_intersect(vector<vec3> &pts, vector<int> &fidxs, vec3 &o, vec3 &d) {
    if(bvh == NULL) { cerr << "ray_intersect: no BVH" << endl; return; } // no BVH, return false
    // Find the triangles who's BV's intersect the ray.  
    vector<mesh_triangle *> res;  
    bvh->Query(res, o.x(), o.y(), o.z(), d.x(), d.y(), d.z());
    if(res.size() == 0) return; 
    // See if the ray actually intersects the triangle.
    pts.reserve(res.size());
    fidxs.reserve(res.size());
    for(size_t i = 0; i < res.size(); i++) {
      vec3 v; 
      if(res[i]->ray_intersect(v, o, d)) {
	pts.push_back(v);
	fidxs.push_back(res[i]->fidx);
      }
    }
  }

  // Returns true on ray intersection. n contains the intersection
  // point and fidx has the face index of the intersected triangle.
  bool mesh::ray_intersect(vec3 &n, int &fidx, vec3 &o, vec3 &d) {
    if(bvh == NULL) return false; // no BVH, return false

    vector<vec3> pts; vector<int> fidxs;
    ray_intersect(pts, fidxs, o, d);
    if(pts.size() == 0) return false;
 
    // Find closest intersection point on interseted triangles relative
    // to the given origin.
    float mindsq = distance3sq(pts[0], o);  
    n = pts[0];
    fidx = fidxs[0];
    for(size_t p = 1; p < pts.size(); p++) {
      float dsq = distance3sq(o, pts[p]);
      if(dsq < mindsq) { 
	// Update return values with closest intersection point.
	n = pts[p]; fidx = fidxs[p]; 
	mindsq = dsq; 
      }
    }
    return true;
  }

  // NOTE: voxelize() assumes that no vertex is outside of the
  // dimension of the volume.  NOTE: Voxelize should return a volume
  // similar to the corresponding funtion in volume8::triangulate().
  void mesh::voxelize(volume8 &v, uint8 bg, uint8 fg, bool min_dim) {
    v.fill(bg); // Clear the volume.

    // Get triangles of the mesh surface.
    vector<triangle> triangles; as_triangles(triangles);

    if(min_dim) {
      // Shifts triangles by minimum dimension.
      vec3 v0(x_min, y_min, z_min); 
      for(size_t t = 0; t < triangles.size(); t++) {
	triangles[t].v[0] -= v0;
	triangles[t].v[1] -= v0;
	triangles[t].v[2] -= v0;
      }
    }

    vector<ivec3> res;
    for(size_t t = 0; t < triangles.size(); t++) {
      res.clear(); v.tri_intersect(res, triangles[t]);
      for(size_t r = 0; r < res.size(); r++) v(res[r]) = fg;
    }

    // TODO: This functionality can get a little tricky. The topology
    // can get complex and we need to consider inside/outside
    // relations. Also there are challenges with self-intersecting
    // meshes!!!  We want to AND those properly as in computational
    // solid geometry (CSG) ops.  For now I do what's simplest and AND
    // self intersecting mesh regions.

    // Find components of dark bg regions in voxel data.
    vector< vector<ivec3> > comps; v.components(comps, 0, 0);
    for(size_t c = 0; c < comps.size(); c++) {
      vector<ivec3> &lc = comps[c];
      // Do not fill component with foreground value if is against the
      // edge of the volume.
      bool fill_component = true;
      for(size_t i = 0; i < lc.size(); i++) {
	if(lc[i].x() == 0 || lc[i].x() == v.width - 1  ||
	   lc[i].y() == 0 || lc[i].y() == v.height - 1 ||
	   lc[i].z() == 0 || lc[i].z() == v.depth - 1) {
	  fill_component = false; break;
	}
      }
      if(fill_component) { for(size_t i = 0; i < lc.size(); i++) v(lc[i]) = 255; }
    }
  }

  void mesh::sample(vector<uint8> &v_samples, volume8 &v) {
    // Get triangles for this mesh.
    vector<triangle> triangles; as_triangles(triangles);
    v_samples.resize(triangles.size());
    for(size_t i = 0; i < v_samples.size(); i++) v_samples[i] = 0;
    vector<ivec3> res;
    for(size_t fidx = 0; fidx < triangles.size(); fidx++) {
      // Intersect triangle with cubes. 
      res.clear(); v.tri_intersect(res, triangles[fidx]);
      // Compute average value of intersected voxels.
      float avg = 0, N = res.size();
      for(size_t r = 0; r < res.size(); r++) {
	avg += (float)v(res[r].x(), res[r].y(), res[r].z()) / N;
      }
      v_samples[fidx] = avg;
    }
  }

  void mesh::surface_bandpass(vector<int> &bandpass, vector<uint8> &v_samples, 
			      int fg, int bg, int T, int blurIter) {
    vector<uint8> v2 = v_samples;
    blur_samples(v2, blurIter);

    // Uses convention that facelabel is set to 1 if non-edge region,
    // 0 if an edge region.
    bandpass.resize(v2.size());
    for(size_t i = 0; i < v2.size(); i++) {
      if((int)v_samples[i] - (int)v2[i] > T) bandpass[i] = fg;
      else bandpass[i] = bg;
    }
    //vector<int> tmp = facelabels;
    //dilate_label(tmp, 3, 0, 1);
    //facelabels = tmp;
  }
  
  void mesh::blur_samples(vector<uint8> &src, int maxIter) {
    vector<uint8> tmp(src.size());
    for(int iter = 0; iter < maxIter; iter++) {
      // Compute the average value of the current face and its
      // neighbors given the input sample image.
      for(size_t fidx = 0; fidx < faces.size(); fidx++) {
	mesh_face &f = faces[fidx];
	float avg = (float)src[fidx] / 4.0f;
	for(int j = 0; j < 3; j++) avg += (float)src[f.adj[j]] / 4.0f;
	tmp[fidx] = avg;
      }
      src = tmp;
    }
  }

  void mesh::dilate_label(vector<int> &dst, vector<int> &src, int niter, int fg, int bg) {
    dst = src;
    vector<int> tmp(dst.size());
    for(int iter = 0; iter < niter; iter++) {
      for(size_t fidx = 0; fidx < faces.size(); fidx++) {
	if(dst[fidx] == fg) tmp[fidx] = fg; // fg face, leave alone
	else if(dst[fidx] == bg) {
	  // If a bg face, count the number of neighboring fg faces.
	  mesh_face &f = faces[fidx];
	  int cnt = 0;
	  for(int j = 0; j < 3; j++) { if(f.adj[j] >= 0 && dst[f.adj[j]] == fg) cnt++; }
	  if(cnt > 0) tmp[fidx] = fg; else tmp[fidx] = bg;
	}
      }
      dst = tmp;
    }
  }

  void mesh::dilate_fg(vector<int> &dst, vector<int> &src, int niter, int fg) {
    dst = src;
    vector<int> tmp(dst.size());
    for(int iter = 0; iter < niter; iter++) {
      // Scan through faces.
      tmp = dst; 
      for(size_t fidx = 0; fidx < faces.size(); fidx++) {
	if(dst[fidx] == fg) {
	  mesh_face &f = faces[fidx];
	  // Label adjacent faces with fg value.
	  for(int j = 0; j < 3; j++) { if(f.adj[j] >= 0) tmp[f.adj[j]] = fg;  }
	}
      }
      dst = tmp;
    }
  }


  void mesh::fill_bg(vector<int> &dst, vector<int> &src, int niter, int bg) {
    dst = src;
    int sel_idx = 0;
    vector<int> tmp(dst.size());
    for(int iter = 0; iter < niter; iter++) {
      tmp = dst; 
      for(size_t fidx = 0; fidx < faces.size(); fidx++) {
	if(dst[fidx] != bg) continue;
	// If a bg face, count the number of neighboring fg faces.
	mesh_face &f = faces[fidx];
	int fg_cnt = 0;
	for(int j = 0; j < 3; j++) { if(f.adj[j] >= 0 && dst[f.adj[j]] != bg) fg_cnt++; }
	if(fg_cnt == 1) {
	  // Single fg. Update current with value.
	  for(int j = 0; j < 3; j++) {
	    if(f.adj[j] >= 0 && dst[f.adj[j]] != bg) { tmp[fidx] = dst[f.adj[j]]; break; } 
	  }
	}
	else if(fg_cnt >= 2) {  // >=2 fg faces, pick max count
	  if(f.adj[0] >= 0  && f.adj[1] >= 0 && dst[f.adj[0]] == dst[f.adj[1]]) {
	    tmp[fidx] = dst[f.adj[0]];
	  }
	  else if(f.adj[0] >= 0  && f.adj[2] >= 0 && dst[f.adj[0]] == dst[f.adj[2]] ) {
	    tmp[fidx] = dst[f.adj[0]];
	  }
	  else if(f.adj[1] >= 0  && f.adj[2] >= 0 && dst[f.adj[1]] == dst[f.adj[2]] ){
	    tmp[fidx] = dst[f.adj[1]];	      
	  }
	  else {
	    // Otherwise pick one arbitrarily (cycle selection). Don't have to worry about 
	    // inf looping here since fg_cnt >= 2 assures no f.adj[sel_idx] == -1.
	    while(true) {
	      if(f.adj[sel_idx] < 0) { sel_idx++; if(sel_idx == 3) sel_idx = 0; }
	      else {
		tmp[fidx] = dst.at(f.adj[sel_idx]);
		sel_idx++; if(sel_idx == 3) sel_idx = 0;
		break;
	      }
	    }
	    tmp[fidx] = dst[f.adj[0]];	      
	  }
	}
      }
      dst = tmp;
    }
  }

  // Watershed segment voxel samples.  I is an image with intensity
  // information.  O is coded as follows: 0 = background label, <= -1
  // leave alone label, >= 1, label to grow by Watershed
  void mesh::watershed(vector<uint8> &I, vector<int> &O) {
    map<int, queue<int> > PQ;

    for(int p = 0; p < num_faces(); p++) {
      if(O[p] > 0) { // > 0 label to use
	// Put all marker/seed pixels with background neighbors into
	// the priority queue.
	bool bgNeighbor = false;
	mesh_face &f = faces[p];      
	for(int i = 0; i < 3; i++) {
	  if(f.adj[i] < 0) continue;
	  int q = f.adj[i];
	  if(O[q] == 0) bgNeighbor = true; //  0 is the back ground label
	}
	if(bgNeighbor) PQ[(int)I[p]].push(p); // Separate queues for each intensity.
      }
    }

    // NOTE: 0 == unlabeled
    while(!PQ.empty()) {
      // Get the queue for the lowest intensity level.
      int vimg = PQ.begin()->first; 
      queue<int> curQ = PQ.begin()->second;
      PQ.erase( PQ.begin() );
      while(curQ.empty() == false) { // Empty this queue.
	int p = curQ.front(); curQ.pop();
	// Propagate the label of p to all unlabeled neighbors of p.
	mesh_face &f = faces[p];
	for(int i = 0; i < 3; i++) {
	  if(f.adj[i] < 0) continue;
	  int q = f.adj[i];
	  if(O[q] == 0) { // Check to see if q is an unlabeled neighbor.
	    O[q] = O[p]; // Propagate the label of p to unlabeled neighbors.
	    // Is the intensity at the propagated position less than
	    // (or equal to) the current queue intensity? Add to
	    // current queue. This deals "reasonably" with plateaus in
	    // the image.
	    if( I[q] <= vimg ) curQ.push(q);
	    else PQ[I[q]].push(q); // Otherwise add to corresponding intensity queue.
	  }
	}
      }
    }
  }

  void mesh::boundary_edges(vector<int> &facelabels, vector<vec3> &v0, vector<vec3> &v1) {
    v0.clear(); v1.clear();
    for(size_t fidx = 0; fidx < faces.size(); fidx++) {
      mesh_face &f = faces[fidx];
      for(int a = 0; a < 3; a++) {
	if(f.adj[a] >= 0 && facelabels[fidx] != facelabels[f.adj[a]]) {
	  int i = -1, j = -1;
	  switch(a) {
	  case 0: i = 0; j = 1; break; // edge 0 = vert 0 and 1
	  case 1: i = 0; j = 2; break; // edge 1 = vert 0 and 2
	  case 2: i = 1; j = 2; break; // edge 2 = vert 1 and 2
	  }
	  v0.push_back(vertices[f.vert[i]].pos);
	  v1.push_back(vertices[f.vert[j]].pos);
	}
      }
    }
  }
  void mesh::boundary_edges(vector<vec3> &v0, vector<vec3> &v1) {
    v0.clear(); v1.clear();
    for(size_t fidx = 0; fidx < faces.size(); fidx++) {
      mesh_face &f = faces[fidx];
      for(int a = 0; a < 3; a++) {
	if(f.adj[a] < 0) {
	  int i = -1, j = -1;
	  switch(a) {
	  case 0: i = 0; j = 1; break; // edge 0 = vert 0 and 1
	  case 1: i = 0; j = 2; break; // edge 1 = vert 0 and 2
	  case 2: i = 1; j = 2; break; // edge 2 = vert 1 and 2
	  }
	  v0.push_back(vertices[f.vert[i]].pos);
	  v1.push_back(vertices[f.vert[j]].pos);
	}
      }
    }
  }

  // Dijkstra's algorithm. Single source shortest path. 
  void mesh::Dijkstra(vector<int> &previous, vector<float> &dist, int source) {
    const float inf = numeric_limits<float>().max();
    dist.resize(num_faces());
    previous.resize(num_faces());

    util::PQi<float> Q(num_faces(), dist);
    for(int fidx = 0; fidx < num_faces(); fidx++) {
      dist[fidx] = inf;
      previous[fidx] = -1;
    }
    dist[source] = 0;
    // Fill the priority queue.
    for(int fidx = 0; fidx < num_faces(); fidx++) Q.insert(fidx);

    while(!Q.empty()) {
      int u = Q.extractMin();
      if(dist.at(u) == inf) break; // No more accessible faces.

      // For each adjacency.
      for(int i = 0; i < 3; i++) {
	int v = faces[u].adj[i]; if(v < 0) continue;
	float alt = dist[u] + distance3sq(face_centroid(u), face_centroid(v));
	if(alt < dist[v]) {
	  dist[v] = alt; Q.decreaseKey(v);
	  previous[v] = u;
	}
      }
    }
  }
  void mesh::ShortestPath(vector<int> &path, vector<int> &previous, int target) {
    int u = target;
    // Path will end up with triangle centroids.
    path.clear();
    while(previous[u] >= 0) {
      path.push_back(u);
      u = previous[u];
    }
    path.push_back(u);
    reverse(path.begin(), path.end());
  }

  void mesh::edge_vert(vector<vec3> &edge_vert) {
    edge_vert.clear();
    set<int> edge_vidxs;
    for(int fidx = 0; fidx < num_faces(); fidx++) {
      if(faces[fidx].adj[0] < 0) {
	edge_vidxs.insert(faces[fidx].vert[0]);
	edge_vidxs.insert(faces[fidx].vert[1]);
      }
      if(faces[fidx].adj[1] < 0) {
	edge_vidxs.insert(faces[fidx].vert[0]);
	edge_vidxs.insert(faces[fidx].vert[2]);
      }
      if(faces[fidx].adj[2] < 0) {
	edge_vidxs.insert(faces[fidx].vert[1]);
	edge_vidxs.insert(faces[fidx].vert[2]);
      }
    }
    for(set<int>::iterator it = edge_vidxs.begin(); it != edge_vidxs.end(); it++) {
      edge_vert.push_back(vertices.at(*it).pos);
    }
  }

  void mesh::pca_dims(vector<float> &dims, float alpha) {
    // Get mesh verticies, subtract centroid
    vector<vec3> mesh_vert(vertices.size());
    for(size_t v = 0; v < vertices.size(); v++) mesh_vert[v] = vertices[v].pos - centroid;

    // Compute PCA using mesh vertices.
    vector<vec3> pca3; vector<float> pca3_x;
    PCA3(pca3, pca3_x, mesh_vert);

    dims.reserve(3); dims.resize(3);
    for(int i = 0; i < 3; i++) {
      vec3 dir1 = pca3.at(i).max_project(mesh_vert); // Make sure to rescale to micron axes.
      vec3 pca3_neg = -pca3[i];
      vec3 dir2 = pca3_neg.max_project(mesh_vert);
      dir1 *= alpha;
      dir2 *= alpha;
      dims[i] = dir1.length() + dir2.length();
    }
  }

  void mesh::pca_longdir(vec3 &pos1, vec3 &pos2, float len_frac) {
    // Get mesh verticies, subtract centroid
    vector<vec3> mesh_vert(vertices.size());
    for(size_t v = 0; v < vertices.size(); v++) mesh_vert[v] = vertices[v].pos - centroid;

    // Compute PCA using mesh vertices.
    vector<vec3> pca3; vector<float> pca3_x;
    PCA3(pca3, pca3_x, mesh_vert);

    vec3 dir1 = pca3.at(2).max_project(mesh_vert); // Make sure to rescale to micron axes.
    vec3 pca3_neg = -pca3[2];
    vec3 dir2 = pca3_neg.max_project(mesh_vert);

    pos1 = dir1 * len_frac + centroid;
    pos2 = dir2 * len_frac + centroid;

    if(pos1.z() > pos2.z()) swap(pos1, pos2);
  }

  float mesh::bend_measure() {
    if(vertices.size() == 0) return -1;
    // Get mesh verticies, subtract centroid
    vector<vec3> mesh_vert(vertices.size());
    for(size_t v = 0; v < vertices.size(); v++) mesh_vert[v] = vertices[v].pos - centroid;

    // Compute PCA using mesh vertices.
    vector<vec3> pca3; vector<float> pca3_x;
    PCA3(pca3, pca3_x, mesh_vert);
    vec3 pc2 = pca3[2], pc1 = pca3[1];

    vector<float> Rsq_vals(mesh_vert.size());
    for(size_t v = 0; v < mesh_vert.size(); v++) {
      float s2 = pc2.scalar_proj(mesh_vert[v]);
      float s1 = pc1.scalar_proj(mesh_vert[v]);
      vec3 on_plane = s2 * pc2 + s1 * pc1;
      Rsq_vals[v] = geom::distance3sq(on_plane, mesh_vert[v]);
    }

    sort(Rsq_vals.begin(), Rsq_vals.end());
  
    return Rsq_vals.at(0.75 * Rsq_vals.size());
  }

  float mesh::compute_sa(float alpha) {
    float surface_area = 0;
    vector<triangle> tri; as_triangles(tri);    
    for(size_t fidx = 0; fidx < tri.size(); fidx++) {
      triangle t = tri[fidx];
      // Map from isotropic voxel space to micron space.
      t[0] *= alpha; t[1] *= alpha; t[2] *= alpha;
      surface_area += t.area();
    }
    return surface_area;
  }

};

