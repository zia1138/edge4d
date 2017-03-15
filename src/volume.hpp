#ifndef __VOLUME_HPP__
#define __VOLUME_HPP__

#include <vector> 
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <limits>


// TODO: Replace all uses of vec3 with QVector.
#include <QMatrix4x4>

#include "geometry.hpp"

namespace vol {
  using namespace std;
  using namespace geom;

  template< typename T > 
  struct cmp_size_gt { bool operator () (const vector<T> &a, const vector<T> &b) const { return a.size() > b.size();} };

  template< typename T > 
  struct cmp_size_lt { bool operator () (const vector<T> &a, const vector<T> &b) const { return a.size() < b.size();} };

  // Builds a sphereical structural element with radius r.
  inline void sphere_strel(vector<ivec3> &strel, int r) {
    int rsq = r*r;
    for(int z = -r; z <= r; z++) 
      for(int y = -r; y <= r; y++) 
	for(int x = -r; x <= r; x++) 
	  if(x * x + y * y + z * z <= rsq) strel.push_back(ivec3(x,y,z));
  }
  
  typedef unsigned char uint8;   // 8-bit data
  typedef unsigned short uint16; // 16-bit image data

  // Adapted from http://graphics.stanford.edu/~seander/bithacks.html
  // Returns the floor form of binary logarithm for a 32 bit integer
  inline int floorLog2(unsigned int v) {  unsigned r = 0;  while (v >>= 1) r++; return r; }
  // Returns true of integer is a lower of 2. 
  inline bool isPow2(unsigned int v) { return v && !(v & (v - 1)); }

  // Converts pixel/voxel radius to sigma for Gaussian blur (adapted
  // from G.I.M.P.)  NOTE: This formula is related to the fullwidth
  // maximum heigh of a Guassin.  Basically it for a pixel radius set
  // at the the < 1% (or 1/255) boundary of a gaussin this gives you
  // the corresponding sigma.
  inline float px2sigma(float px) { return sqrt(-(px * px / (2 * log(1.0/255.0)))); }

  // Intensity histogram with bins. 
  template <typename T = uint16> struct hist_t {
    vector<long> hist;  // counts of values in histogram bins 
    vector<double> cdf;  // cumulative density function 
    T vmin, vmax, range;
    // Initialize thistogram.
    hist_t() { 
      hist.resize(numeric_limits<T>::max()+1, 0);  
      cdf.resize(numeric_limits<T>::max()+1,1);  
      vmin = vmax = range = 0; 
    }
    // Add a value to the histogram. 
    void add(T v) { hist[v]++; }
    // Compute cumulative density function (CDF). 
    void compute_cdf() {
      double N = 0, cum = 0;
      for(size_t i = 0; i < hist.size(); i++) N += hist[i];
      cout << "N = " << N << endl;
      size_t b = 0; 
      bool show_med = false;
      while(b < hist.size()) { 
	cum += hist[b]; 
	cdf[b] = double(cum) / double(N); 
	if(cdf[b] > 0.5 && !show_med) { cout << "median = " << cdf[b] << " " << b  << endl; show_med = true; }
	b++; 
      }
    }
    // This histogram adjustment assumes the image values are in the range 0 to 1.
    void adjust(double fmin, double fmax) {
      compute_cdf();
      // Compute new max and min values and range based on bin position.
      vmin = 0; 
      while(vmin < numeric_limits<T>::max() && cdf[vmin] < fmin) vmin++;
      vmax = numeric_limits<T>::max(); 
      while(vmax > 0 && 1.0f - cdf[vmax] < fmax) vmax--;
      cout << (long)vmin << " rr " << (long)vmax << endl;
      if(vmax <= vmin) { cerr << "hist_t:adjust(), vmax <= vmin" << endl; exit(1); }
      range = vmax - vmin;
    }
    // Adjust histogram values.
    float map(T v) { 
      if(v > vmax) return numeric_limits<T>::max(); else if(v < vmin) return 0;  
      else return (float(v) - float(vmin)) / float(range) * (float)numeric_limits<T>::max(); 
    }
  };

  // Stores a volume of histograms of an original image stack.
  template <typename T> struct histvol_t {
    int width, height, depth;
    float fmin, fmax; // Desired min and max frequency for adjustment.
    hist_t<T> global;    // global histogram (of entire volume)
    vector< hist_t<T> > data; // data stores all of the histograms
    histvol_t(int w, int h, int d, float fmin_, float fmax_) { 
      fmin = fmin_; fmax = fmax_; data.resize(w*h*d); 
      width = w; height = h; depth = d;
    }
    hist_t<T> & v(int x, int y, int z) { 
      return data.at(z * (width * height) + y * width + x);  
    }
    hist_t<T> & operator () (int x, int y, int z) { return v(x,y,z); }

    // Trilinear interpolation of histograms adjusted values.  This prevents aliasing.
    float map(float xd, float yd, float zd, float val) {
      if(width == 1 && height == 1 && depth == 1) return data[0].map(val);
      int x = (int)xd, y = (int)yd, z = (int)zd;
      int xp1 = x+1, yp1 = y+1, zp1 = z+1;
      float dx = xd - x, dy = yd - y, dz = zd - z;

      // There is no interpolation performed at x = width - 1, y = height - 1, and z = depth-1.
      if(xp1 >= width)  xp1 = x; if(yp1 >= height) yp1 = y; if(zp1 >= depth)  zp1 = z;

      float sample = 
	v(x,y,z).map(val) * (1-dx) * (1-dy) * (1-dz) + // 0 0 0
	v(xp1,y,z).map(val) * dx * (1-dy) * (1-dz) +   // 1 0 0
	v(x,yp1,z).map(val) * (1-dx) * dy * (1-dz) +   // 0 1 0
	v(x,y,zp1).map(val) * (1-dx) * (1-dy) * dz +   // 1 1 0
	v(xp1,y,zp1).map(val) * dx * (1-dy) * dz +     // 1 0 1
	v(x,yp1,zp1).map(val) * (1-dx) * dy * dz +     // 0 1 1
	v(xp1,yp1,z).map(val) * dx * dy * (1-dz) +     // 0 0 1
	v(xp1,yp1,zp1).map(val) * dx * dy * dz;        // 1 1 1

      return sample;
    }
    // Pre-process data in histogram volumes.
    void pre_process() { for(size_t i = 0; i < data.size(); i++) data[i].adjust(fmin,fmax); }
  };

  template <class T> struct volumeT {
    int width, height, depth, odim;
    T *data;
    volumeT() { data = NULL; width = height = depth = 0; }
    volumeT(int w_, int h_, int d_) {
      width = w_; height = h_; depth = d_;  
      data = new T[width * height * depth]; 
      setodim(width, height, depth); 
    }
    virtual ~volumeT() { if(data != NULL) delete[] data; }
    volumeT( volumeT<T> &v2 ) {
      width = v2.width; height = v2.height; depth = v2.depth;
      data = new T[width * height * depth];       
      setodim(width, height, depth);
      copy(v2);
    }
    // Create a copy of the image volume.
    void copy(volumeT<T> &src) {
      long whd = width * height * depth;
      for(long i = 0; i < whd; i++) data[i] = src.data[i];
      setodim(width, height, depth); // update octtree dim
    }
    void pad(int r, T bg) {
      volumeT<T> p(width + 2 * r, height + 2 * r, depth + 2 * r);
      p.fill(bg);
      for(int z = 0; z < depth; z++) 
	for(int y = 0; y < height; y++) 
	  for(int x = 0; x < width; x++) p(x+r,y+r,z+r) = v(x,y,z);	  
      // Copy padded block. 
      delete data; data = p.data; p.data = NULL; 
      // Update dimensions to included pad region.
      width = p.width; height = p.height; depth = p.depth;
    }      

    void unpad(int r) {
      volumeT<T> p(width - 2 * r, height - 2 * r, depth - 2 * r);
      for(int z = 0; z < p.depth; z++) 
	for(int y = 0; y < p.height; y++) 
	  for(int x = 0; x < p.width; x++) p(x,y,z) = v(x+r,y+r,z+r);	  
      // Copy padded block. 
      delete data; data = p.data; p.data = NULL; 
      // Update dimensions to included pad region.
      width = p.width; height = p.height; depth = p.depth;
    }

    void pad_AB(int r, T bg) {
      volumeT<T> p(width, height, depth + 2 * r);
      p.fill(bg);
      for(int z = 0; z < depth; z++) 
	for(int y = 0; y < height; y++) 
	  for(int x = 0; x < width; x++) p(x,y,z+r) = v(x,y,z);	  
      // Copy padded block. 
      delete data; data = p.data; p.data = NULL; 
      // Update dimensions to included pad region.
      width = p.width; height = p.height; depth = p.depth;
    }
    void unpad_AB(int r) {
      volumeT<T> p(width, height, depth - 2 * r);
      for(int z = 0; z < p.depth; z++) 
	for(int y = 0; y < p.height; y++) 
	  for(int x = 0; x < p.width; x++) p(x,y,z) = v(x,y,z+r);	  
      // Copy padded block. 
      delete data; data = p.data; p.data = NULL; 
      // Update dimensions to included pad region.
      width = p.width; height = p.height; depth = p.depth;
    }

    // odim is set to the surrounding octree structure.
    void setodim(int w_, int h_, int d_) {
      unsigned int w = w_, h = h_, d = d_;
      // Make the width, height and depth a power of 2 for octtree recursion.
      if(isPow2(w)) w = 1 << floorLog2(w); else w = 1 << (floorLog2(w) + 1);
      if(isPow2(h)) h = 1 << floorLog2(h); else h = 1 << (floorLog2(h) + 1);
      if(isPow2(d)) d = 1 << floorLog2(d); else d = 1 << (floorLog2(d) + 1);
      odim = max(max(w,h),d);
    }
    // Accessor functions.
    T &v(int x, int y, int z) { return data[z * (width * height) + y * width + x]; }
    T &v(const ivec3 &p) { return v(p[0], p[1], p[2]); }

    // Adapted from http://www.paulinternet.nl/?page=bicubic
    float cubicInterpolate (float p[4], float x) {
      return p[1] + 0.5f * x*(p[2] - p[0] + x*(2.0f*p[0] - 5.0f*p[1] + 4.0f*p[2] - p[3] + x*(3.0f*(p[1] - p[2]) + p[3] - p[0])));
    }
    float bicubicInterpolate (float p[4][4], float x, float y) {
      float arr[4];
      arr[0] = cubicInterpolate(p[0], y); arr[1] = cubicInterpolate(p[1], y); 
      arr[2] = cubicInterpolate(p[2], y); arr[3] = cubicInterpolate(p[3], y);
      return cubicInterpolate(arr, x);
    }
    float tricubicInterpolate (float p[4][4][4], float x, float y, float z) {
      float arr[4];
      arr[0] = bicubicInterpolate(p[0], y, z); arr[1] = bicubicInterpolate(p[1], y, z);
      arr[2] = bicubicInterpolate(p[2], y, z); arr[3] = bicubicInterpolate(p[3], y, z);
      return cubicInterpolate(arr, x);
    }
    T tricubic(float xd, float yd, float zd) {
      int x = (int)xd, y = (int)yd, z = (int)zd;
      float dx = xd - x, dy = yd - y, dz = zd - z;
      // Is this right? 
      float p[4][4][4];
      if(x <= 0 || y <= 0 || z <= 0 || x >= width-3 || y >= height-3 || z >= depth-3) {
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j < 4; j++)
	    for (int k = 0; k < 4; k++) {
	      int px = x + i - 1, py = y + j - 1, pz = z + k - 1;
	      if(px < 0) px = 0; if(py < 0) py = 0; if(pz < 0) pz = 0;
	      if(px >= width) px = width-1; if(py >= height) py = height-1;
	      if(pz >= depth) pz = depth-1;
	      p[i][j][k] = v(px, py, pz);
	    }
      }
      else {
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j < 4; j++)
	    for (int k = 0; k < 4; k++) p[i][j][k] = v(x + i - 1, y + j - 1, z + k - 1);
      }
      float I = tricubicInterpolate(p, dx, dy, dz);
      // Clip interoplated value.
      const T minT = numeric_limits<T>::min(), maxT = numeric_limits<T>::max();
      T IT = I; if(I < (float)minT) IT = minT; if (I > (float)maxT) IT = maxT;
      return IT;
    }
    // Nearest neighbor interpolation.
    T nearest(float xd, float yd, float zd) {
      int x = roundf(xd), y = roundf(yd), z = roundf(zd);
      // Catch out of bounds case, and return value.
      if(x >= width) x = width - 1; if(y >= height) y = height - 1; if(z >= depth) z = depth - 1;
      if(x < 0) x = 0; if(y < 0) y = 0; if(z < 0) z = 0;
      return v(x,y,z);
    }
    T & operator () (int x, int y, int z = 0) { return v(x,y,z); } 
    T & operator () (const ivec3 &p) { return v(p[0], p[1], p[2]); }
    void fill(T s) { for(int i = 0; i < width * height * depth; i++) data[i] = s; }

    // Scale image using tricubic interpolation.
    void scale_tricubic(volumeT<T> &src, int zstart = 0, int zstep = 1) {
      float sx = float(src.width) / float(width), sy = float(src.height) / float(height), sz = float(src.depth) / float(depth); 
      for(int vz = zstart; vz < depth; vz += zstep) 
	for(int vy = 0; vy < height; vy++) 
	  for(int vx = 0; vx < width; vx++) v(vx,vy,vz) = src.tricubic(vx * sx, vy * sy, vz * sz);
    }
    void scale_nearest(volumeT<T> &src, int zstart = 0, int zstep = 1) {
      float sx = float(src.width) / float(width), sy = float(src.height) / float(height), sz = float(src.depth) / float(depth); 
      for(int vz = zstart; vz < depth; vz += zstep) 
 	for(int vy = 0; vy < height; vy++) 
 	  for(int vx = 0; vx < width; vx++) v(vx,vy,vz) = src.nearest(vx * sx, vy * sy, vz * sz);
    }

    void mirrorx(volumeT<T> &src) {
      for(int z = 0; z < depth; z++)
	for(int y = 0; y < height; y++)
	  for(int x = 0; x < width; x++) {
	    v(x,y,z) = src((width - 1) - x, y, z);
	  }
    }


    void rotate90(volumeT<T> &src) {
      for(int z = 0; z < depth; z++)
	for(int y = 0; y < height; y++)
	  for(int x = 0; x < width; x++) {
	    int zold = (width - 1) - x;
	    if(zold >= depth || zold < 0) continue;
	    int xold = z;
	    if(xold >= width || xold < 0) continue;	    
	    v(x,y,z) = src(xold, y, zold);
	  }
    }

    void rotate90minus(volumeT<T> &src) {
      for(int z = 0; z < depth; z++)
	for(int y = 0; y < height; y++)
	  for(int x = 0; x < width; x++) {
	    int zold = x;
	    if(zold >= depth || zold < 0) continue;
	    int xold = (depth - 1) - z;
	    if(xold >= width || xold < 0) continue;	    
	    v(x,y,z) = src(xold, y, zold);
	  }
    }

    double correlation(volumeT<T> &src) {
      double corr = 0;
      for(int z = 0; z < depth; z++)
	for(int y = 0; y < height; y++)
	  for(int x = 0; x < width; x++)  corr += v(x,y,z) * src(x,y,z);
      return corr;
    }

    void ridged(volumeT<T> &src, QMatrix4x4 &transform) {
      for(int z = 0; z < depth; z++)
	for(int y = 0; y < height; y++)
	  for(int x = 0; x < width; x++) {
	    QVector4D p = transform.map(QVector4D(x,y,z, 1));
	    v(x,y,z) = src.nearest(p.x(), p.y(), p.z());
	  }
    }
    
    // x stays the same, y and z update.
    //v(x,y,z) = src(, y, (width-1) - x);
    //v(x,y,z) = src((width - 1) - x, y, z);
    
    // Use octtree recursion to limit voxels interseced with vertex
    // ray.  (x,y,z) and (x+dim,y+dim,z+dim) define octtree cube at
    // current recursion.
    void oct_isect(vector<ivec3> &res, vec3 &p, vec3 &ni, float t, int dim, int x, int y, int z) {
      if(x >= width || y >= height || z >= depth) return; // Out of volume (is this right or is above right?)
      AABB cube(x,x+dim, y,y+dim, z,z+dim);
      bool isect = cube.contains(p[0],p[1],p[2]);
      if(!isect) {
	float lmin, lmax;
	isect = cube.ray_intersect(lmin, lmax, p[0],p[1],p[2], ni[0],ni[1],ni[2]);
	if(isect && fabsf(lmin) > t) isect = false;
      }
      if(isect) {
	if(dim == 1) res.push_back(ivec3(x,y,z));
	else {
	  // 8-way octree recursion
	  int d = dim / 2;
	  oct_isect(res, p, ni, t, d, x    , y    , z    );  // 0 0 0
	  oct_isect(res, p, ni, t, d, x    , y    , z + d);  // 0 0 1
	  oct_isect(res, p, ni, t, d, x    , y + d, z    );  // 0 1 0
	  oct_isect(res, p, ni, t, d, x    , y + d, z + d);  // 0 1 1
	  oct_isect(res, p, ni, t, d, x + d, y    , z    );  // 1 0 0
	  oct_isect(res, p, ni, t, d, x + d, y    , z + d);  // 1 0 1
	  oct_isect(res, p, ni, t, d, x + d, y + d, z    );  // 1 1 0
	  oct_isect(res, p, ni, t, d, x + d, y + d, z + d);  // 1 1 1
	}
      }
    }
    // Returns voxel positions that intersect a given ray at position
    // p with direction n. Voxels must be less than given distance.
    void ray_intersect(vector<ivec3> &res, vec3 p, vec3 n, float dist) { 
      vec3 ninv(1.0f/n[0],1.0f/n[1],1.0f/n[2]); oct_isect(res, p, ninv, dist, odim, 0, 0, 0); 
    }
    // This octtree recursion accelerates ray-voxel intersections. 
    void oct_isect(vector<ivec3> &res, triangle &tri, AABB &tribox, int dim, int x, int y, int z) {
      if(x >= width || y >= height || z >= depth) return; // Out of volume (is this right or is above right?)
      AABB cube(x,x+dim, y,y+dim, z,z+dim);
      if(cube.tri_intersect(tri, tribox)) {
	if(dim == 1) res.push_back(ivec3(x,y,z));
	else {
	  // 8-way octree recursion
	  int d = dim / 2;
	  oct_isect(res, tri, tribox, d, x    , y    , z    );  // 0 0 0
	  oct_isect(res, tri, tribox, d, x    , y    , z + d);  // 0 0 1
	  oct_isect(res, tri, tribox, d, x    , y + d, z    );  // 0 1 0
	  oct_isect(res, tri, tribox, d, x    , y + d, z + d);  // 0 1 1
	  oct_isect(res, tri, tribox, d, x + d, y    , z    );  // 1 0 0
	  oct_isect(res, tri, tribox, d, x + d, y    , z + d);  // 1 0 1
	  oct_isect(res, tri, tribox, d, x + d, y + d, z    );  // 1 1 0
	  oct_isect(res, tri, tribox, d, x + d, y + d, z + d);  // 1 1 1
	}
      }
    }    
    void tri_intersect(vector<ivec3> &res, triangle &tri) {
      // Precompute bounding box of triangle.
      AABB tribox(tri.min_x(), tri.max_x(), tri.min_y(), tri.max_y(), tri.min_z(), tri.max_z());
      oct_isect(res, tri, tribox, odim, 0, 0, 0); 
    }
    // Subtract a volume (of the same size) from the current volume. 
    void subtract(volumeT<T> &v2) {
      long wdh = width * height * depth;
      for(long i = 0; i < wdh; i++) data[i] -= v2.data[i];
    }

    void nhood6_unsafe(vector<ivec3> &hood, ivec3 p) {
      hood.clear();
      ivec3 adj[6];
      adj[0] = ivec3(p[0] - 1, p[1], p[2]); adj[1] = ivec3(p[0] + 1, p[1], p[2]);
      adj[2] = ivec3(p[0], p[1] - 1, p[2]); adj[3] = ivec3(p[0], p[1] + 1, p[2]);
      adj[4] = ivec3(p[0], p[1], p[2] - 1); adj[5] = ivec3(p[0], p[1], p[2] + 1);
      for(int i = 0; i < 6; i++) hood.push_back(adj[i]);
    }
    // Return a neighborhood of voxel positions.
    void nhood6(vector<ivec3> &hood, ivec3 p) {
      hood.clear();
      ivec3 adj[6];
      adj[0] = ivec3(p[0] - 1, p[1], p[2]); adj[1] = ivec3(p[0] + 1, p[1], p[2]);
      adj[2] = ivec3(p[0], p[1] - 1, p[2]); adj[3] = ivec3(p[0], p[1] + 1, p[2]);
      adj[4] = ivec3(p[0], p[1], p[2] - 1); adj[5] = ivec3(p[0], p[1], p[2] + 1);
      for(int i = 0; i < 6; i++) {
	ivec3 &n = adj[i];
	bool outside = n[0] < 0      || n[1] < 0       || n[2] < 0       || 
	  n[0] >= width || n[1] >= height || n[2] >= depth;
	if(!outside) hood.push_back(n);
      }
    }
    void nhood26_unsafe(vector<ivec3> &hood, ivec3 p) {
      hood.clear();
      for(int dx = -1; dx <= 1; dx++) 
	for(int dy = -1; dy <= 1; dy++)
	  for(int dz = -1; dz <= 1; dz++) {
	    if(dx != 0 || dy != 0 || dz != 0) {
	      hood.push_back(ivec3(p[0] + dx, p[1] + dy, p[2] + dz));
	    }
	  }

    }
    void nhood26(vector<ivec3> &hood, ivec3 p) {
      hood.clear();
      for(int dx = -1; dx <= 1; dx++) 
	for(int dy = -1; dy <= 1; dy++)
	  for(int dz = -1; dz <= 1; dz++) {
	    if(dx != 0 || dy != 0 || dz != 0) {
	      ivec3 n(p[0] + dx, p[1] + dy, p[2] + dz);
	      bool outside = n[0] < 0      || n[1] < 0       || n[2] < 0       || 
		n[0] >= width || n[1] >= height || n[2] >= depth;
	      if(!outside) hood.push_back(n);
	    }
	  }
    }
    // Fill histogram with intensity values.    
    void PopulateHistVol(histvol_t<T> &hv) {
      double sx = hv.width / double(width), sy = hv.height / double(height), sz = hv.depth / double(depth);
      for(int z = 0; z < depth; z++) 
	for(int y = 0; y < height; y++) 
	  for(int x = 0; x < width; x++) {
	    // populate global & local histograms
	    hv.global.add(v(x,y,z)); 
	    hv(x * sx, y * sy, z * sz).add(v(x,y,z));
	  }
    }
    void ApplyHistogram(histvol_t<T> &hv, int ystart = 0, int ystep = 1) {
      double sx = hv.width / double(width), sy = hv.height / double(height), sz = hv.depth / double(depth);
      for(int z = 0; z < depth; z++) 
	for(int y = ystart; y < height; y += ystep) 
	  for(int x = 0; x < width; x++) 
	    v(x,y,z) = hv.map(x * sx, y * sy, z * sz, v(x,y,z));
    }

  };


  // 16-bit (single channel) image volume. Primarily used for loading images 
  struct volume16 : public volumeT<uint16> {
    bool is8bitOriginal;
    volume16() { is8bitOriginal = false; }
    volume16(int w, int h, int d) : volumeT<uint16> (w,h,d) { }
    // Load volume from a hyperstack.
    bool load(string filename, int slices, int num_channels, int frame, int channel, int z_depth = 0, int z_start = 0);
  };

  enum thintype_t { thinLINE, thinSURFACE };
  enum thindir_t { thinU = 0, thinD = 1, thinN = 2, thinS = 3, thinE = 4, thinW = 5 }; 

  // Most of the image processing is done using single channel 8-bit image volumes.
  struct volume8 : public volumeT<uint8> { 
    int x0, y0, z0; // Position of the binary volume in a "larger" binary volume.
    volume8() { data = NULL; x0 = y0 = z0 = 0; }
    volume8(int w, int h, int d = 1) : volumeT<uint8>(w,h,d) { x0 = y0 = z0 = 0; fill(0); }
    // Simple threshold, set to FG.
    volume8(volume8 &v, int threshold) : volumeT<uint8>(v.width, v.height, v.depth) {
      x0 = y0 = z0 = 0;
      for(int i = 0; i < (width * height * depth); i++) data[i] = v.data[i] > threshold ? 255 : 0;
    }
    volume8(volume16 &v) : volumeT<uint8>(v.width, v.height, v.depth) {
      x0 = y0 = z0 = 0;
      if(v.is8bitOriginal) { for(int i = 0; i < (width * height * depth); i++) data[i] = v.data[i]; }
      else {
	for(int i = 0; i < (width * height * depth); i++) data[i] = (float(v.data[i]) / 65535.0f) * 255.0f;
      }
    }
    // Copy constructor.
    volume8(volume8 &v) : volumeT<uint8>(v.width, v.height, v.depth) {
      x0 = y0 = z0 = 0;
      for(long i = 0; i < (width * height * depth); i++) data[i] = v.data[i];
    }
    // Creates a binary volume using given foreground voxel positions.
    volume8(vector<ivec3> &fg, int border = 2) {
      border = max(border, 1); // border must be 1 or greater
      data = NULL; x0 = y0 = z0 = 0; // Remember the "lower left bottom" corner.
      if(fg.size() > 0) {
	// Get minimum and maximum (x,y,z) voxel positions.
	ivec3 minpos, maxpos; 
	geom::bounding_box(minpos, maxpos, fg);
	x0 = minpos.x(); y0 = minpos.y(); z0 = minpos.z();
	// Add 2 x border to width and height.
	// NOTE: Added a +1 here.. make sure it doesn't affect anything!!!!
	// NOTE: Added a +1 here.. make sure it doesn't affect anything!!!!
	// NOTE: Added a +1 here.. make sure it doesn't affect anything!!!!
	// NOTE: Added a +1 here.. make sure it doesn't affect anything!!!!
	width =  (maxpos.x() - minpos.x() + 1) + 2 * border; 
	height = (maxpos.y() - minpos.y() + 1) + 2 * border; 
	depth =  (maxpos.z() - minpos.z() + 1) + 2 * border;
	data = new uint8[width * height * depth]; 
	fill(0); // Zero the volume
	// Update pixels within volume.
	for(size_t i = 0; i < fg.size(); i++) {
	  // Shift by border.
	  int x = (fg[i].x() - x0) + border, y = (fg[i].y() - y0) + border, z = (fg[i].z() - z0) + border;
	  v(x,y,z) = 255;
	}
	// Adjust upper left top corner position by border.
	x0 -= border; y0 -= border; z0 -= border;
      }
      else {
	cerr << "volume8: fg.size() == 0" << endl;
	x0 = y0 = z0 = -border;
	width = 2*border; height = 2*border; depth = 2*border;
	data = new uint8[width*height*depth];
	fill(0);
      }
    }

    // Find connected components by depth first search.
    void dfs(volumeT<char> &labels, ivec3 start, vector<ivec3> &comp, uint8 fg = 255, bool use_nhood26 = false);
    void components(vector< vector<ivec3> > &comps, int minvxls, uint8 fg, bool use_nhood26 = false);
    // Applies shift to resulting component.
    void components2(vector< vector<ivec3> > &comps, int minvxls = 0, uint8 fg = 255, bool use_nhood26 = false) {
      components(comps, minvxls, fg, use_nhood26);
      ivec3 shift0(x0,y0,z0);
      for(size_t c = 0; c < comps.size(); c++) 
	for(size_t i = 0; i < comps[c].size(); i++)  comps[c][i] += shift0;
    } 

    // Finds voxels neighboring a background voxel (outlining) fg component.
    void outline_vox(vector<ivec3> &ovox, uint8 fg = 255, bool use_nhood26 = false);

    // Uses Marching Cubes to triangulate a binary volume.  See
    // critical details in code to assure geom::mesh() returns the same volume.
    void triangulate(vector<vec3> &vtable, vector<face> &ftable, 
		     int mcdim, uint8 fg = 255, uint8 bg = 0);

    void dilate(vector<ivec3> &update, int r, uint8 fg, uint8 bg, 
		float alpha = 0.0, bool invert = false, int zstart = 0, int zstep = 1);
    void dilate(int r, uint8 fg, uint8 bg, float alpha = 0.0) {
      vector<ivec3> upd;
      dilate(upd, r, fg, bg, alpha);
      for(size_t i = 0; i < upd.size(); i++) v(upd[i]) = fg;
    }
    void substrel(vector<int> &hist, vector<ivec3> &s, int x, int y, int z);
    void addstrel(vector<int> &hist, vector<ivec3> &s, int x, int y, int z);

    void set(vector<ivec3> &update, uint8 val) { for(size_t i = 0; i < update.size(); i++) v(update[i]) = val; }
  
    // Inverts the image volume. 
    void invert(uint8 fg = 255, uint8 bg = 0) {
      for(int i = 0; i < width * height * depth; i++) {
	if(data[i] == fg) data[i] = bg;	else if(data[i] == bg) data[i] = fg;
      }
    }

    void subtract_and_clip(uint8 subval) {
      for(int z = 0; z < depth; z++) 
	for(int y = 0; y < height; y++) 
	  for(int x = 0; x < width; x++) v(x,y,z) = max((int)v(x,y,z) - (int)subval, 0);
    }

    void invert_grayscale() {
      for(int i = 0; i < width * height * depth; i++) data[i] = 255 - data[i];
    }

    // Clears a border along the edge of the volume.
    void clear_border(int brdr, uint8 val = 0) {
      for(int z = 0; z < depth; z++)
	for(int y = 0; y < height; y++)
	  for(int x = 0; x < width; x++) {
	    if(x < brdr || y < brdr || z < brdr ||
	       x >= width - brdr || y >= height - brdr || z >= depth - brdr)
	      v(x,y,z) = val;
	  }
    }

    // Clears a border along XY dimensions only.
    void clear_borderXY(int brdr, uint8 val = 0) {
      for(int z = 0; z < depth; z++)
	for(int y = 0; y < height; y++)
	  for(int x = 0; x < width; x++) {
	    if(x < brdr || y < brdr || x >= width - brdr || y >= height - brdr) v(x,y,z) = val;
	  }
    }


    void subtract_and_clip(volume8 &v1, volume8 &v2, int threshold) {
      long wdh = width * height * depth;
      for(long i = 0; i < wdh; i++) {
	int delta = (int)v1.data[i] - (int)v2.data[i];
	data[i] = delta > threshold ? 255 : 0;
      }
    }

    // O(rwdh) algorithm for rank filtering, using stacked histograms.
    void rankfilt(volume8 &dst, int r, float alpha = 0.5, int zstart = 0, int zstep = 1);
    // Helper functions for median computation.
    inline void addhist(int hist[256], int med, int &cum, vector<ivec3> &s, int x, int y, int z);
    inline void subhist(int hist[256], int med, int &cum, vector<ivec3> &s, int x, int y, int z);

    // "A Sequential 3D Thinning Algorithm and Its Medical Applications" Kalman Palagyi et. al
    bool is_border_point(thindir_t direction, int x, int y, int z);
    bool is_endpoint(int x, int y, int z, uint8 fg = 255);
    bool is_surface_endpoint(int x, int y, int z, uint8 fg = 255);
    bool is_finished(thintype_t thintype, int x, int y, int z, uint8 fg = 255);
    void thinptl(vector<ivec3> &ptl, thindir_t direction, thintype_t thintype, int zstart = 0, int zstep = 1);
    int thin_apply(vector<ivec3> &ptl, thintype_t thintype);
    int thin_subiter(thindir_t direction, thintype_t thintype);
    void thin_init(); // Clears one voxel border.
    // Single threaded driver call of thinning.
    void thin(thintype_t thintype = thinLINE);

    // Fills holes with fg value, all but the largest background region.
    void fill_holes(uint8 fg = 255, uint8 bg = 0);
    
    // Gaussain blur in X, Y, and Z dimensions.
    void GaussX(float sigmaX, int ystart = 0, int ystep = 1); 
    void GaussY(float sigmaY, int xstart = 0, int xstep = 1); 
    void GaussZ(float sigmaZ, int xstart = 0, int xstep = 1);

    // Computes the feature transform using the algorithm in
    // "A Linear Time Algorithm for Computing Exact Euclidean Distance
    // Transforms of Binary Images in Arbitrary Dimensions" IEEE
    // Transactions on Pattern Analysis and Machine Intelligence. Vol
    // 25, No 2., Februrary 2003, p265-270.
    // The feature transform gives the (x,y,z) position of the nearest
    // foreground voxel in Euclidean distance.
    int sqr(int x) { return x * x; }
    bool RemoveFT(int d2u, int d2v, int d2w, int ud, int vd, int wd);
    void VoronoiFT(volumeT<ivec3> &F, int d, ivec3 j);
    void ComputeFT(volumeT<ivec3> &F, uint8 fg, int d, ivec3 j);
    void ComputeFT(volumeT<ivec3> &F, uint8 fg) { ComputeFT(F, fg, 3, ivec3()); }


    void below_threshold(volume8 &bv, int t, uint8 bg = 0, uint8 fg = 255) {
      long whd = width * height * depth;
      for(long i = 0; i < whd; i++) bv.data[i] = data[i] <= t ? fg : bg;
    }

    // Save 8-bit volume as a .TIFF file.
    bool save(string filename);
    bool save_mha(string filename, float dx, float dy, float dz); 
  };

  // Algorithms from Vincent, L., "Morphological Grayscale Reconstruction in Image Analysis: Applications and Efficient
  // Algorithms," IEEE Transactions on Image Processing, Vol. 2, No. 2, April, 1993, pp. 176-201.
  void binary_reconstruct(volume8 &J, volume8 &I, uint8 fg = 255, bool use_nhood26 = false);
  void reconstruct_raster(volume8 &J, volume8 &I, bool use_nhood26 = false);

  struct Hmatch_t {
    float dist;  
    ivec3 shift; 
    int max_label, max_overlap;
    bool operator < (const Hmatch_t &b) const { return dist < b.dist; }
  };

  // 32-bit (signed) integer volume, used for Euclidean Distance Transform and Watershed segmentation.
  struct volume32 : public volumeT<int> {
    volume32() : volumeT<int>() { }
    volume32(int w, int h, int d = 1) : volumeT<int> (w,h,d) { fill(0); }

    // Accessor functions for EDT.
    int &D(int x, int y, int z) { return v(x,y,z); }
    int &D(const ivec3 &p) { return v(p[0], p[1], p[2]); }
    // Single thread EDT computation.
    void ComputeEDT(volume8 &I, uint8 fg = 255) { ComputeEDT(I, fg, 3, ivec3()); }
    // Compute EDT for x and y dimensions (useful for multithreading).
    void ComputeEDTstep1(volume8 &I, uint8 fg, int zstart = 0, int zstep = 1) {
      for(int z = zstart; z < depth;  z+=zstep) ComputeEDT(I, fg, 2, ivec3(-1,-1,z) );
    }
    // Finally compute it for z dimension (useful for multithreading).
    void ComputeEDTstep2(int xstart = 0, int xstep = 1) {
      for(int x = xstart; x < width; x+=xstep)
	for(int y = 0; y < height; y++) VoronoiEDT(3, ivec3(x,y,-1));
    }

    // "A Linear Time Algorithm for Computing Exact Euclidean Distance
    // Transforms of Binary Images in Arbitrary Dimensions" IEEE
    // Transactions on Pattern Analysis and Machine Intelligence. Vol
    // 25, No 2., Februrary 2003, p265-270.
    int sqr(int x) { return x * x; }
    bool RemoveEDT(int d2u, int d2v, int d2w, int ud, int vd, int wd);
    void VoronoiEDT(int d, ivec3 j);
    void ComputeEDT(volume8 &I, uint8 fg, int d, ivec3 j);

    void above_threshold(volume8 &bv, int t, uint8 bg = 0, uint8 fg = 255) {
      long whd = width * height * depth;
      for(long i = 0; i < whd; i++) bv.data[i] = data[i] >= t ? fg : bg;
    }
    void below_threshold(volume8 &bv, int t, uint8 bg = 0, uint8 fg = 255) {
      long whd = width * height * depth;
      for(long i = 0; i < whd; i++) bv.data[i] = data[i] <= t ? fg : bg;
    }
    void hystersis_threshold(volume8 &hout, int thresh_min, int thresh_max);

    // Seeds regions, starting with a label of 1 ...
    void seed(vector< vector<ivec3> > &seeds);

    // Negates all values in 32-bit signed volume.
    void negate() { long whd = width * height * depth; for(long i = 0; i < whd; i++) data[i] = -data[i]; }

    // If the given volume has foreground value, set the data value to -1.
    void seed(volume8 &bv, uint8 fg, int val) {
      long whd = width * height * depth;
      for(long i = 0; i < whd; i++) if(bv.data[i] == fg) data[i] = val;
    }
    // Finds components with value greater than zero. 
    void components(vector< vector<ivec3> > &final, bool compact_empty = true);

    // Accesor using a O() and an ivec3.
    int &O(const ivec3 &p) { return v(p[0], p[1], p[2]); }

    // Seeded watershed algorithm. Propagates labeled regions (which
    // must be have value > 0) using given image (usually negated
    // EDT).  Use negative values in I to exclude regions uses 0 to
    // designate unlabeld regions.
    void watershed(volume32 &I, bool use_nhood26 = false);

    // Fractional hausdorf distance of given fg voxel set. Note
    // volume32 must contain the EDT. Also returns -1 if distance
    // can't be computed.  voxel positions can be outside volume.
    float frac_hausdorff(const vector<ivec3> &vox, float alpha = 0.8);

    // Template matching using fractional Hausdorff.
    float frac_hausdorff_match(ivec3 &min_shift, vector<ivec3> &vox, int d, float alpha = 0.8);

    // Get value of region with maximum overlap.
    void max_overlap(int &max_label, int &max_cnt, vector<ivec3> &vox);

    // Returns best hasudorff matches with unique maximum overlap
    // regions, according to the given labels. Results sorted for
    // minimum distance to max.
    void frac_hausdorff_match(vector<Hmatch_t> &matches, 
			      vector<ivec3> &vox, vector<ivec3> &voxfill, 
			      int d, float alpha, volume32 &labels);

    // Save 32-bit volume in a RAW format.
    void save_raw(const string &filename);

    bool load_mha_float(const string &filename);
  };

  void grayscale_reconstruct(volume32 &J, volume32 &I, bool use_nhood26 = false);

  // Computes shortest path in image volume.
  float ShortestPath(vector<ivec3> &path, volume8 &v, ivec3 src, ivec3 targ, uint8 fg = 255);

  // Generates a synthetic data set.
  void GenSynth(volume8 &v, int num_cells, float max_dist);

  // Current version of segmentation algorithm.
  void HSegment3_Nucs(volume8 &edge_vol, volume8 &v, 
		      vector< vector<ivec3> > &nuc_vox, 
		      int dmin, int dmax, 
		      int internal_cnt,  
		      vector< vector<ivec3> > &final, 
		      int min_comp, int max_comp, int noise_comp);

  ivec3 hausdorff_align(vector<ivec3> &A, vector<ivec3> &B);

}

#endif

