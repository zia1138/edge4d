#ifndef __MESH__HPP__
#define __MESH__HPP__

#include <vector>
#include <iostream>
#include <limits>

// define is needed for VBO function calls (makes triangle drawing a lot faster).
#define GL_GLEXT_PROTOTYPES

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glext.h>
#else
#include <GL/gl.h>
#include <GL/glext.h>
#endif

#include "geometry.hpp"
#include "volume.hpp"

namespace geom {
  using namespace vol;

  struct mesh_vertex { vec3 pos, normal; };
  struct mesh_face {
    mesh_face() { adj[0] = adj[1] = adj[2] = -1; vert[0] = vert[1] = vert[2] = -1; }
    vec3 normal; 
    int adj[3];  // edge 0 = vertex 0 and 1, edge 1 = vertex 0 and 2, edge 2 = vertex 1 and 2
    int vert[3]; // indices (in the vertex array) of all vertices (mesh_vertex)
  };

  // A triangle with a face index. Used in the BVH for accelerated intersection calculations.
  struct mesh_triangle : public triangle { int fidx; };

  struct mesh {
    vector<mesh_face> faces; 
    vector<mesh_vertex> vertices;
    vector< vector<int> > facelist; // adjacent face list for each vertex
    // NOTE: If any of the vertices change position the BVH needs to be rebuilt!!!
    vector<mesh_triangle> tri_face; // mapping from a triangle to a face index
    BVH<mesh_triangle> *bvh;        // BVH on triangles

    vec3 centroid;
    vector<vec3> pca3, axes3; vector<float> pca3_x;

    vector< triangle > tribuf;      // triangle buffer for quick drawing

    // Mesh follows a model where data is copied to VBO after updates.
    bool vboInitFlag; // set to true if VBO is initialized
    GLuint vboId; // vertex buffer with triangles

    int num_faces() { return faces.size(); } int num_vert() { return vertices.size(); }

    vec3 &vertex(int vidx) { return vertices[vidx].pos; }
    
    // Interface and geometry for the bounding volume heirarchy.
    float x_max, y_max, z_max, x_min, y_min, z_min, x_center, y_center, z_center;
    float x() { return x_center; } float y() { return y_center; } float z() { return z_center; }
    float min_x() { return x_min; } float max_x() { return x_max; }
    float min_y() { return y_min; } float max_y() { return y_max; }
    float min_z() { return z_min; } float max_z() { return z_max; }

    void bounding_box(vec3 &mn, vec3 &mx) {
      mn[0] = x_min; mn[1] = y_min; mn[2] = z_min;
      mx[0] = x_max; mx[1] = y_max; mx[2] = z_max;
    }
    
    // Recomputes min_x,y,z max_x,y,z and center_x,y,z. Call this when vertex positions change.
    float compute_sa(float alpha);
    void reset_range(); 
    void compute_axes();
    void build_bvh();
    void rebuild_bvh() { if(bvh != NULL) { delete bvh; bvh = NULL; build_bvh(); } }

    vec3 face_centroid(int fidx) {
      mesh_face &face = faces[fidx];
      vec3 v[3];
      for(int i = 0; i < 3; i++) v[i] = vertices[face.vert[i]].pos;
      return (v[0] + v[1] + v[2]) / 3.0;
    }
    void Dijkstra(vector<int> &previous, vector<float> &dist, int source);
    void ShortestPath(vector<int> &path, vector<int> &previous, int target);
    // Transform face index path into a vector of face centroids.
    void CentroidPath(vector<vec3> &path, vector<int> &fidx_path) {
      path.clear();
      for(size_t i = 0; i < fidx_path.size(); i++) 
	path.push_back(face_centroid(fidx_path[i]));
    }

    void edge_faces(vector<int> &edge_fidx) {
      for(int fidx = 0; fidx < num_faces(); fidx++) {
	if(faces[fidx].adj[0] < 0 || faces[fidx].adj[1] < 0 ||  faces[fidx].adj[2] < 0 ) {
	  edge_fidx.push_back(fidx);
	}
      }
    }
    void edge_vert(vector<vec3> &edge_vert);
    
    // Intersect mesh with ray return nearest intersected face index
    // and intersection point n.
    bool ray_intersect(vec3 &n, int &fidx, vec3 &o, vec3 &d);
    // Returns all intersection points and intersected faces.
    void ray_intersect(vector<vec3> &pts, vector<int> &fidxs, vec3 &o, vec3 &d);

    // Create a mesh from a vertex table and a face table.
    void face_adjacencies(); // F-F adjacencies
    void rebuild_adj();
    void rebuild(vector<vec3> &vtable, vector<face> &ftable);

    // Compute two flavors of normals.
    void compute_face_normals();  void compute_vertex_normals();

    // Find connected components in a mesh.
    void connected_components(vector< vector<int> > &components, int min_tri = 1000); 
    void dfs_label(vector<int> &labels, int start, int label);

    // Given a list of faces returns a sub-mesh as a vertex and face table.
    void as_mesh(vector<int> &subfaces, vector<vec3> &vtable, vector<face> &ftable);

    double compute_volume(double vxd);

    // Return centroid of the mesh.
    vec3 compute_centroid() {
      vec3 C; float N = num_vert();
      for(int vidx = 0; vidx < num_vert(); vidx++) C += vertices[vidx].pos / N;
      return C;
    }

    // Taubin lambda/mu mesh smoothing. Algorithm adapted from TriMesh2.
    void smooth(int niters);
    vec3 smooth_disp(int vidx);

    void smooth_iter(vector<vec3> &disp, float stepsize);

    void volpts(vector<ivec3> &vs);
    // Gets voxels along thin boundary with thickness dthresh.
    void volpts(vector<ivec3> &pts, int dthresh);
    // Gets voxels above and below a given plane.
    void volpts(vector<ivec3> &above, vector<ivec3> &below, plane3 &plane);

    void volpts_internal(vector<ivec3> &pts, ivec3 &min_v);

    // Update mesh meta information, call when vertex positions change.
    void update_meta(bool update_bvh = true) { 
      reset_range();                // change min/max x,y,z ranges
      compute_face_normals();       // update both face normals
      compute_vertex_normals();     // and vertex normals
      if(update_bvh) rebuild_bvh(); // rebuild the BVH since points will have moved

      as_triangles(tribuf);         // update triangle buffer (todo, replace with strips?)
    }

    // Updates GL VBO buffers. NOTE: This can't be called in a thread!!!!
    void update_gl() { // ** NOT THREAD SAFE!!!! **
      // VBO shit here.
      if(vboInitFlag == false) { vboInitFlag = true; glGenBuffersARB(1, &vboId); }
      else {
	glDeleteBuffersARB(1, &vboId); cout << "deleting buffer" << endl;
	glGenBuffersARB(1, &vboId);
      }
      vector<vec3> vert, norm;  
      // Copy data to  vertex buffer. Note that each triangle vertex must have a normal.
      // TODO: Does a graphics card treat these as a vertex normal, interpolating between?
      as_vnorm(vert, norm); 
      glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboId);
      glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(vec3) * (vert.size() + norm.size()), 0, GL_STATIC_DRAW_ARB);
      glBufferSubDataARB(GL_ARRAY_BUFFER_ARB, 0, sizeof(vec3) * vert.size(), &vert[0]);
      glBufferSubDataARB(GL_ARRAY_BUFFER_ARB, sizeof(vec3) * vert.size(), sizeof(vec3) * norm.size(), &norm[0]);
      glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
    }

    // Return all faces as triangles.
    void as_triangles(vector<triangle> &tri) {
      tri.resize(faces.size());
      for(int fidx = 0; fidx < (int)tri.size(); fidx++) tri[fidx] = as_triangle(fidx); 
    }
    // Get a sub-set of triangles.  
    void as_triangles(vector<int> &subfaces, vector<triangle> &tri) {
      tri.resize(subfaces.size());
      for(size_t i = 0; i < tri.size(); i++) tri[i] = as_triangle(subfaces[i]);
    }

    void as_vnorm(vector<vec3> &vert, vector<vec3> &norm) {
      vert.resize(3*faces.size());
      norm.resize(3*faces.size());
      for(size_t fidx = 0; fidx < faces.size(); fidx++) {
	mesh_face &face = faces[fidx];
	int fidx3 = 3*fidx;
	vert[fidx3]   = vertices[face.vert[2]].pos;
	vert[fidx3+1] = vertices[face.vert[1]].pos;
	vert[fidx3+2] = vertices[face.vert[0]].pos;
	norm[fidx3]   = face.normal; norm[fidx3+1] = face.normal; norm[fidx3+2] = face.normal;
      }
    }
    
    // Returns a given face index as a triangle.
    triangle as_triangle(int fidx) {
      triangle t;
      mesh_face &face = faces[fidx];
      t.normal = face.normal;
      t.v[0] = vertices[face.vert[0]].pos; t.v[1] = vertices[face.vert[1]].pos; t.v[2] = vertices[face.vert[2]].pos;
      return t;
    }

    mesh(vector<vec3> &vtable, vector<face> &ftable) { 
      vboId = 0; vboInitFlag = false;
      bvh = NULL; 
      rebuild(vtable, ftable); 
    }
    ~mesh() { 
      if(bvh != NULL) delete bvh; 
      if(vboInitFlag) glDeleteBuffersARB(1, &vboId);
    }

    // Converts mesh into a binary volume. 
    void voxelize(volume8 &v, uint8 bg = 0, uint8 fg = 255, bool min_dim = false);

    // Move vertex positions to match image volume.
    vec3 refine_disp(volume8 &v, int vidx, float dist, float sign = -1);
    void refine(volume8 &v, float dist, int max_iter, float alpha, float stepsize);

    // Sample voxels into image volume (converts 8-bit to 16-bit samples).
    void sample(vector<uint8> &vsamples, volume8 &v);

    // Operations for image processing on surface labels.
    void surface_threshold();
    void surface_bandpass(vector<int> &bandpass, vector<uint8> &v_samples, 
			  int fg = 1, int bg = 0, int T = 100, int blurIter = 120);
    void blur_samples(vector<uint8> &src, int maxIter = 120);
    void dilate_label(vector<int> &dst, vector<int> &src, 
		      int niter = 3, int fg = 1, int bg = 0);

    void fill_bg(vector<int> &dst, vector<int> &src, int niter, int bg);
    void dilate_fg(vector<int> &dst, vector<int> &src, int niter, int fg);

    void clean_small_comps(vector<int> &facelabels, int num_comps, int thresh = 30);
    void dfs_paint(vector<int> &facelabels, int start, int label, int dthresh = 1);
    int label_components(vector<int> &labels);
    // Performs marker assisted watershed segmentation on the mesh.
    void watershed(vector<uint8> &I, vector<int> &O);

    // Finds vertices of bounday edges for a set of face labels.
    void boundary_edges(vector<int> &facelabels, vector<vec3> &v0, vector<vec3> &v1);

    // Finds edges of triangles that aren't adjacent to another triangle.
    void boundary_edges(vector<vec3> &v0, vector<vec3> &v1);

    void component_areas(vector<int> &facelabels, float alpha, vector<float> &areas);

    void pca_dims(vector<float> &dims, float alpha);
    void pca_longdir(vec3 &pos1, vec3 &pos2, float len_frac);

    float bend_measure();

  };

};

#endif // __MESH__HPP__
