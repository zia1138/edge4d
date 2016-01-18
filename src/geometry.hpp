#ifndef __GEOMETRY_HPP__
#define __GEOMETRY_HPP__

#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>

#include <math.h>

// Moller code in tribox3.cpp TODO: Get rid of this and use a simpler
// (but slower) approach to triangle-AABB intersections.
int triBoxOverlap(float boxcenter[3],float boxhalfsize[3],float triverts[3][3]);

namespace geom {
  using namespace std;

  // Return determinant of 2x2 matrix 
  inline float det2(float a, float b,
		    float c, float d) { return (a * d - b * c); }

  // Return determinant of 3x3 matrix.
  inline float det3(float a, float b, float c, 
		    float d, float e, float f, 
		    float g, float h, float i) {
    return (a * (e * i - h * f) - b * (d * i - g * f) + c * (d * h - g * e));
  }

  // A simple 3x3 matrix.
  struct matrix3 {
    float data[3][3];
    float &v(int i, int j) { return data[i][j]; }
    float v(int i, int j) const { return data[i][j]; }
    matrix3() { 
      v(0,0) = 0; v(0,1) = 0; v(0,2) = 0;
      v(1,0) = 0; v(1,1) = 0; v(1,2) = 0;
      v(2,0) = 0; v(2,1) = 0; v(2,2) = 0;
    }
    float &operator () (int i, int j) {  return v(i,j); }
    float operator () (int i, int j) const {  return v(i,j); }
    float trace() { return v(0,0) + v(1,1) + v(2,2); }
    // Compute the determinant.
    float det() const {
      return det3(v(0,0), v(0,1), v(0,2),
		  v(1,0), v(1,1), v(1,2),
		  v(2,0), v(2,1), v(2,2));
    }
    // Multiply two 3x3 matricies.
    matrix3 operator * (const matrix3 &p) {
      matrix3 R;
      for(int i = 0; i < 3; i++)  for(int j = 0; j < 3; j++) R(i,j) = v(i,0) * p(0,j) + v(i,1) * p(1,j) + v(i,2) * p(2,j);
      return R;
    }
    // Matrix subtract. 
    matrix3 operator - (const matrix3 &p) {
      matrix3 R; for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++) R(i,j) = v(i,j) - p(i,j);
      return R;
    }
  };

  //  Left-multiply a 3x3 matrix by a scalar.
  inline matrix3 operator * (float s, const matrix3 &A) {
    matrix3 R;
    R(0,0) = s * A(0,0); R(0,1) = s * A(0,1); R(0,2) = s * A(0,2);
    R(1,0) = s * A(1,0); R(1,1) = s * A(1,1); R(1,2) = s * A(1,2);
    R(2,0) = s * A(2,0); R(2,1) = s * A(2,1); R(2,2) = s * A(2,2);
    return R;
  }

  // Return a 3x3 identity matrix.
  inline matrix3 eye3() {
    matrix3 A; 
    A(0,0) = 1; A(0,1) = 0; A(0,2) = 0;
    A(1,0) = 0; A(1,1) = 1; A(1,2) = 0;
    A(2,0) = 0; A(2,1) = 0; A(2,2) = 1;
    return A;
  }

  // Compute an inverse 3x3 matrix. 
  matrix3 inv(const matrix3 &A); 

  // A simple 3-d vector class.
  struct vec3 {
    float v[3];
    vec3() { v[0] = v[1] = v[2] = 0; }
    vec3(float v0, float v1, float v2) { v[0] = v0; v[1] = v1; v[2] = v2; }
    vec3(float v2[3]) { v[0] = v2[0];  v[1] = v2[1]; v[2] = v2[2]; }
    vec3(const vec3 &z) { v[0] = z[0]; v[1] = z[1]; v[2] = z[2]; }
    vec3(float x) { v[0] = v[1] = v[2] = x; }
    void zero() { v[0] = v[1] = v[2] = 0; }

    void operator = (float x) { v[0] = x; v[1] = x; v[2] = x; }

    float &operator[] (int dim) { return v[dim]; }
    float operator[] (int dim) const { return v[dim]; }
  
    vec3 operator-() { return vec3(-v[0],-v[1],-v[2]);  }

    float &x() { return v[0]; } float &y() { return v[1]; } float &z() { return v[2]; }

    vec3 operator * (float a) { return vec3(v[0] * a, v[1] * a, v[2] * a);  }
    vec3 operator * (const vec3 &z) { return vec3(v[0] * z[0], v[1] * z[1], v[2] * z[2]); }
    vec3 operator + (const vec3 &z) { return vec3(v[0] + z[0], v[1] + z[1], v[2] + z[2]); }
    vec3 operator - (const vec3 &z) { return vec3(v[0] - z[0], v[1] - z[1], v[2] - z[2]); }
    void operator += (const vec3 &z) { v[0] += z[0]; v[1] += z[1]; v[2] += z[2]; }
    void operator *= (const vec3 &z) { v[0] *= z[0]; v[1] *= z[1]; v[2] *= z[2]; }
    void operator += (float a) { v[0] += a; v[1] += a; v[2] += a; }
    void operator -= (float a) { v[0] -= a; v[1] -= a; v[2] -= a; }
    void operator *= (float a) { v[0] *= a; v[1] *= a; v[2] *= a; }
    void operator /= (float a) { v[0] /= a; v[1] /= a; v[2] /= a; }
    void operator -= (const vec3 &z) { v[0] -= z[0]; v[1] -= z[1]; v[2] -= z[2]; }
    vec3 operator / (float a) { return vec3(v[0] / a, v[1] / a, v[2] / a); }

    vec3 operator % (const vec3 &z) {   // Cross product
      vec3 x;
      x[0] = (v[1]*z[2]) - (v[2]*z[1]);
      x[1] = (v[2]*z[0]) - (v[0]*z[2]);
      x[2] = (v[0]*z[1]) - (v[1]*z[0]);
      return x;
    }

    float l2normsq() { return v[0] * v[0] + v[1] * v[1] + v[2] * v[2]; }
    float length() { return sqrt(l2normsq()); }
    float l2norm() { return length(); }
    void normalize() {  float Z = length(); v[0] /= Z; v[1] /= Z; v[2] /= Z;  }
    float dot(const vec3 &z) { return v[0] * z[0] + v[1] * z[1] + v[2] * z[2]; }
    float angle(const vec3 &b) { 
      float dotprod = dot(b);
      if(dotprod > 1 || dotprod < -1) return 0; // Addresses numeric problems.
      else return acosf(dotprod); 
    }
    float abs_angle(const vec3 &b) { 
      float dotprod = fabs(dot(b));
      if(dotprod > 1) return 0; // Addresses numeric problems.
      else return acosf(dotprod); 
    }

    // Magnitude of b in the direction of this vector.
    float scalar_proj(const vec3 &b) { return dot(b) / l2norm(); }

    // Vector projection of b onto this vector. 
    vec3 proj(const vec3 &b) { return vec3(v[0], v[1], v[2]) * (dot(b) / l2normsq()); }

    vec3 max_mag_project(vector<vec3> &vs) {
      if(vs.size() == 0) return vec3(0,0,0);
      float max_scalar = fabsf(scalar_proj(vs[0]));
      vec3 maxv = vs[0];
      for(size_t i = 1; i < vs.size(); i++) {
	float scalar = fabsf(scalar_proj(vs[i]));
	if(scalar > max_scalar) { max_scalar = scalar; maxv = vs[i];}
      }
      return proj(maxv);
    }
    vec3 max_project(vector<vec3> &vs, int *max_idx_ref = NULL) {
      if(vs.size() == 0) return vec3(0,0,0);
      int max_idx = 0;
      float max_scalar = scalar_proj(vs[0]);
      vec3 maxv = vs[0];
      for(size_t i = 1; i < vs.size(); i++) {
	float scalar = scalar_proj(vs[i]);
	if(scalar > max_scalar) { max_scalar = scalar; maxv = vs[i]; max_idx = i; }
      }
      if(max_idx_ref != NULL) *max_idx_ref = max_idx;
      return proj(maxv);
    }
  };

  inline float rad2deg(float rad) { return rad * (180.0 / M_PI); }

  // Integer 3d vector class.  TODO: Maybe make vec3 template typed
  // <T> to minimize code repeated.
  struct ivec3 {
    ivec3() { v[0] = v[1] = v[2] = 0; }
    ivec3(int x, int y, int z) { v[0] = x; v[1] = y; v[2] = z; }
    ivec3(const ivec3 &z) { v[0] = z[0]; v[1] = z[1]; v[2] = z[2]; }    
    int v[3];
    int &x() { return v[0]; } int &y() { return v[1]; } int &z() { return v[2]; }
    int &operator[] (int dim) { return v[dim]; }
    int operator[] (int dim) const { return v[dim]; }
    bool operator == (const ivec3 &z) const { return v[0] == z[0] && v[1] == z[1] && v[2] == z[2]; }
    bool operator != (const ivec3 &z) const { return v[0] != z[0] || v[1] != z[1] || v[2] != z[2]; }
    ivec3 operator + (const ivec3 &z) { return ivec3(v[0] + z[0], v[1] + z[1], v[2] + z[2]); }
    ivec3 operator - (const ivec3 &z) { return ivec3(v[0] - z[0], v[1] - z[1], v[2] - z[2]); }
    ivec3 operator + (int s) { return ivec3(v[0] + s, v[1] + s, v[2] + s); }
    ivec3 operator - (int s) { return ivec3(v[0] - s, v[1] - s, v[2] - s); }

    void operator += (int s) { v[0] += s; v[1] += s; v[2] += s; }
    void operator -= (int s) { v[0] -= s; v[1] -= s; v[2] -= s; }

    void operator += (const ivec3 &z) { v[0] += z[0]; v[1] += z[1]; v[2] += z[2]; }
    void operator -= (const ivec3 &z) { v[0] -= z[0]; v[1] -= z[1]; v[2] -= z[2]; }

    void operator = (int x) { v[0] = x; v[1] = x; v[2] = x; }

    vec3 to_vec3() { vec3 p; p[0] = v[0]; p[1] = v[1]; p[2] = v[2]; return p; }
    int l2normsq() { return v[0] * v[0] + v[1] * v[1] + v[2] * v[2]; }
  };

  inline void bounding_box(ivec3 &mn, ivec3 &mx, const vector<ivec3> &vs) {
    mn = numeric_limits<int>::max();
    mx = -(numeric_limits<int>::max());
    for(size_t i = 0; i < vs.size(); i++) {
      for(int d = 0; d < 3; d++) mn[d] = min(vs[i][d], mn[d]);
      for(int d = 0; d < 3; d++) mx[d] = max(vs[i][d], mx[d]);
    }
  }
  inline void bounding_box(vec3 &mn, vec3 &mx, const vector<ivec3> &vs) {
    ivec3 mn_i, mx_i;
    bounding_box(mn_i, mx_i, vs);
    mn = mn_i.to_vec3();  mx = mx_i.to_vec3();
  }


  inline float l2norm(vec3 v) { return v.l2norm(); }

  // Compute a 3x3 matrix as an outer product of two vectors.
  inline matrix3 outer3(vec3 &a, vec3 &b);

  // Left-multiply a vector by a scalar.
  inline vec3 operator * (float a, const vec3 &v) { return vec3(a * v[0] , a * v[1], a * v[2]); }
  inline vec3 operator * (const matrix3 &A, const vec3 &v) {
    vec3 r;
    r[0] = A(0,0) * v[0] + A(0,1) * v[1] + A(0,2) * v[2];
    r[1] = A(1,0) * v[0] + A(1,1) * v[1] + A(1,2) * v[2];
    r[2] = A(2,0) * v[0] + A(2,1) * v[1] + A(2,2) * v[2];
    return r;
  }

  // An axis with an origin and 3 directions.
  struct axis { vec3 origin, d[3]; };

  // Compute PCA of a a set of points.  x has the eigenvalues from
  // smallest to largest.
  vec3 Centroid3(vector<vec3> &points);
  vec3 Centroid3(vector<ivec3> &points);
  bool PCA3(vector<vec3> &v, vector<float> &x, vector<vec3> &points);

  // Clears largest component. If there are multiple components with the max size, all are cleared.
  void ClearLargest(vector< vector<ivec3> > &comps);

  // Compute the distance between two vectors.
  inline float distance3sq(vec3 s, vec3 t) { 
    float dx = s[0] - t[0], dy = s[1] - t[1], dz = s[2] - t[2];
    return dx * dx + dy * dy + dz * dz;
  }
  inline float distance3sq(ivec3 s, ivec3 t) { 
    float dx = s[0] - t[0], dy = s[1] - t[1], dz = s[2] - t[2];
    return dx * dx + dy * dy + dz * dz;
  }
  inline float distancesq_xy(vec3 s, vec3 t) { 
    float dx = s[0] - t[0], dy = s[1] - t[1];
    return dx * dx + dy * dy;
  }


  inline float distance3(vec3 s, vec3 t) { return sqrtf( distance3sq(s,t) ); }

  // Rodrigues' rotation formula 
  inline vec3 rotate(vec3 &v, vec3& a, float theta) { 
    return v * cosf (theta) + a * (v.dot(a)) * (1 - cosf (theta)) - (v % a) * sinf(theta);
  }

  // Compute triangle area in 3-d.
  float area(vec3 &a, vec3 &b, vec3 &c);

  inline vec3 ComputeNormal(vec3 &p0, vec3 &p1, vec3 &p2) {
    vec3 v1 = p1 - p0, v2 = p2 - p0;
    vec3 normal = v1 % v2;
    normal.normalize();
    return normal;
  }

  // An independent 3-d triangle with an intersection test.
  struct triangle { 
    vec3 v[3], normal;
    //triangle() { color[1] = 1.0; } // Green is the default.

    void as_plane(vec3 &n, float &d) {  n = normal;  d = -n.dot(v[0]); }
    float area() { return geom::area(v[0], v[1], v[2]); }
    vec3 centroid() { return v[0] / 3.0 + v[1] / 3.0 + v[2] / 3.0; }

    // Triangle centroid.
    float x() const { return (1.0/3.0) * (v[0][0] + v[1][0] + v[2][0]); } 
    float y() const { return (1.0/3.0) * (v[0][1] + v[1][1] + v[2][1]); } 
    float z() const { return (1.0/3.0) * (v[0][2] + v[1][2] + v[2][2]); } 

    float minf(const float a, const float b) const { return a < b ? a : b; }
    float maxf(const float a, const float b) const { return a > b ? a : b; }

    // Accessors are used by BVH.
    float max_x() const { return maxf(v[0][0], maxf(v[1][0], v[2][0])); }
    float max_y() const { return maxf(v[0][1], maxf(v[1][1], v[2][1])); }
    float max_z() const { return maxf(v[0][2], maxf(v[1][2], v[2][2])); }

    float min_x() const { return minf(v[0][0], minf(v[1][0], v[2][0])); }
    float min_y() const { return minf(v[0][1], minf(v[1][1], v[2][1])); }
    float min_z() const { return minf(v[0][2], minf(v[1][2], v[2][2])); }

    void ComputeNormal() { 
      vec3 v1 = v[1] - v[0], v2 = v[2] - v[0]; 
      normal = v1 % v2; 
      normal.normalize();
    }
    vec3 &operator[] (int tex) { return v[tex]; }
    // "Fast, minimum storage ray-triangle intersection".  Tomas
    // Moller and Ben Trumbore.  Journal of Graphics Tools,
    // 2(1):21--28, 1997.  Adapted from
    // http://www.mathworks.com.au/matlabcentral/fileexchange/25058-raytriangle-intersection
    // by Jesus P. Mena-Chalco
    bool ray_intersect(vec3 &i, vec3 &o, vec3 &d) {
      vec3 e1 = v[1] - v[0], e2 = v[2] - v[0];
      vec3 q = d % e2;
      float a = e1.dot(q);
      if(a > -1e-6f && a < 1e-6f) return false;
      float f = 1.0f/a;
      vec3 s = o - v[0];
      float u = f * s.dot(q);   // u,v barycentric coordinates of intersection
      if(u < 0.0f) return false;
      vec3 r = s % e1;
      float v = f * d.dot(r);   // u,v barycentric coordinates of intersection
      if( v < 0.0f || u + v > 1.0f) return false;
      float t = f * e2.dot(r); // t distance from the ray origin
      i = o + t * d;
      return t > 0.0f; // NOTE: i = o + t * d, intersection = origin + t*direction 
    }

    // TODO: Add triangle triangle intersection.
    bool tri_intersect(triangle &) {
      // TODO: Use ray_intersect above
      // ray 1 = v[0], d = v[1] - v[0], t = distance from ray origin < l2normsq(d)
      // ray 2 = v[1], d = v[2] - v[1], 
      // ray 3 = v[2], d = v[3] - v[2],
      cout << "not implemented yet" << endl;
      return false;
    }
  };

  // A face stores indices of given vertices in a vector<vec3> array.
  struct face { int vert[3]; }; 


  struct plane3 {
    vec3 v; float d;
    // Construct plane from point and normal vector
    void from_pnorm(vec3 p, vec3 n) { n.normalize(); v = n; d = -p.dot(n); }
    float pdist(const vec3 &p) { return v.dot(p) + d; }
    bool neg_side(const vec3 &p) { return pdist(p) < 0; }
    bool pos_side(const vec3 &p) { return pdist(p) > 0; }
  };

  // perspective camera geometry
  struct camera3 {
    vec3 eye; // position
    vec3 towards, right, up; // direction
    float xfov, yfov;  // field of view
    float neardist, fardist; // distance of near clipping plane and far plane

    // ix and iy are integer image positions.
    void create_ray(vec3 &o, vec3 &p, int ix, int iy, int width, int height) {
      float px = ix, py = iy;

      // Place image points in -1 < x < 1 and -1 < y < 1.
      float x = 2.0 * (px + 0.5) / float(width)  - 1.0;
      float y = -(2.0 * (py + 0.5) / float(height) - 1.0); // TODO: Check if flipped!!!!!

      vec3 &T = towards, &U = up, &R = right;

      // Compute direction based on camera geometry.
      /*p = T + U * y * tanf(yfov/2) + R * x * tanf(xfov/2);
	p.normalize();*/

      T.normalize(); U.normalize(); R.normalize();
      float D = 1.0 / tanf(yfov/2); // distance to projection plane
      float ar = float(width) / float(height); // aspect ratio
      p = D * T + U * y + R * ar * x; // ray direction
      p.normalize();
      o = eye + neardist * p; // Move ray origin to the near plane.
    }
    // For all the calls below dx and dy define a delta position in
    // the movement of the mouse. center is the center of the scene.
    void rotate_scene(int dx, int dy, int width, int height, vec3 &center) {
      float vx = (float) dx / (float) width, vy = (float) dy / (float) height;
      float theta = 4.0 * (fabsf(vx) + fabsf(vy));
      vec3 vector = (right * vx) + (up * vy);
      vec3 rotation_axis = towards % vector; rotation_axis.normalize();
      eye -= center; // Don't forget to translate back from center.
      eye = rotate(eye, rotation_axis, theta) + center; 
      towards = rotate(towards, rotation_axis, theta);  up = rotate(up, rotation_axis, theta);
      right = towards % up;
      up = right % towards;
      towards.normalize(); up.normalize(); right.normalize();
    }    
    void scale_scene(int dx, int dy, int width, int height, vec3 &center) {
      float factor = (float) dx / (float) width + (float) dy / (float) height;
      factor = exp(2.0 * factor);
      factor = (factor - 1.0) / factor;
      vec3 translation = (center - eye) * factor; 
      eye += translation;
    }
    void translate_scene(int dx, int dy, int width, int height, vec3 &center) {
      float length = distance3(center, eye) * tanf(yfov);
      float vx = length * (float) dx / (float) width;
      float vy = length * (float) dy / (float) height;
      vec3 translation = -((right * vx) + (up * vy));
      eye += translation;
    }
  };

  // Axis aligned bounding box (AABB) with a ray intersection function.
  struct AABB {
    float x_min, x_max, y_min, y_max, z_min, z_max;

    // Simplest constructor. It just assigns max and min elements as specified.
    AABB(float x_min_, float x_max_, float y_min_, float y_max_, float z_min_, float z_max_){ 
      x_min = x_min_; x_max = x_max_; y_min = y_min_; y_max = y_max_; z_min = z_min_; z_max = z_max_;
    }
    // If no ranges are specified, the region spans the entire (x, y, z) space.
    AABB() { 
      x_min = -(numeric_limits<float>::max()); x_max = numeric_limits<float>::max();
      y_min = -(numeric_limits<float>::max()); y_max = numeric_limits<float>::max();
      z_min = -(numeric_limits<float>::max()); z_max = numeric_limits<float>::max();
    }

    // Adapted from http://ompf.org/ray/ray_box.html Branchless
    // Ray/Box intersections by Mueller/Geimer Also has numeric fixes
    // (for possible NaNs) described on
    // http://www.flipcode.com/archives/SSE_RayBox_Intersection_Test.shtml
    float minf(const float a, const float b) { return a < b ? a : b; }
    float maxf(const float a, const float b) { return a > b ? a : b; }
    // The paramters [ idx, idy, idz ] = [ 1 / dx, 1 / dy , 1 / dz], the inverse of the direction vector.
    // o = [x,y,z] and d = [dx,dy,dz]
    // o + lmin * d = first intersection point
    // o + lmax * d = second intersection point
    bool ray_intersect(float &lmin, float &lmax, float x, float y, float z, float idx, float idy, float idz) { 
      float l1, l2, fl1a, fl1b, fl2a, fl2b;
      l1	= (x_min - x) * idx;
      l2	= (x_max - x) * idx;
      fl1a = minf(l1, numeric_limits<float>::infinity());  fl2a = minf(l2, numeric_limits<float>::infinity());
      fl1b = maxf(l1, -numeric_limits<float>::infinity()); fl2b = maxf(l2, -numeric_limits<float>::infinity());
      float lmin_x = min(fl1b,fl2b);
      float lmax_x = max(fl1a,fl2a);
      
      l1	= (y_min - y) * idy;
      l2	= (y_max - y) * idy;
      fl1a = minf(l1, numeric_limits<float>::infinity());  fl2a = minf(l2, numeric_limits<float>::infinity());
      fl1b = maxf(l1, -numeric_limits<float>::infinity()); fl2b = maxf(l2, -numeric_limits<float>::infinity());
      float lmin_y	= minf(fl1b,fl2b);
      float lmax_y	= maxf(fl1a,fl2a);

      l1	= (z_min - z) * idz;
      l2	= (z_max - z) * idz;
      fl1a = minf(l1, numeric_limits<float>::infinity());  fl2a = minf(l2, numeric_limits<float>::infinity());
      fl1b = maxf(l1, -numeric_limits<float>::infinity()); fl2b = maxf(l2, -numeric_limits<float>::infinity());
      float lmin_z	= minf(fl1b,fl2b);
      float lmax_z	= maxf(fl1a,fl2a);

      lmin = maxf(lmin_z,maxf(lmin_y, lmin_x));
      lmax = minf(lmax_z,minf(lmax_y, lmax_x));

      //return ((lmax > 0.f) && (lmax >= lmin));
      //return ((lmax > 0.f) && (lmax > lmin));
      //return ((lmax >= 0.f) && (lmax >= lmin));
      return ((lmax >= 0.f) && (lmax >= lmin));
    }
    bool ray_intersect(float x, float y, float z, float idx, float idy, float idz) { 
      float lmin, lmax; 
      return ray_intersect(lmin, lmax, x,y,z, idx,idy,idz);
    }
    bool contains(float x, float y, float z) {
      return ( x_min <= x && x <= x_max ) && ( y_min <= y && y <= y_max ) && ( z_min <= z && z <= z_max );
    }

    // Arvo's algorithm from "A Simple Method for Box-Sphere
    // Intersection Testing", by Jim Arvo, in "Graphics Gems",
    // Academic Press, 1990.
    bool sphere_intersect(AABB &bbox, float r, float cx, float cy, float cz) {
      if(intersect_AABB(bbox) == false) return false;  // No bbox intersect, skip.
      float s, d = 0;

      if( cx < x_min ) { s = cx - x_min; d += s * s; } else
      if( cx > x_max ) { s = cx - x_max; d += s * s; }

      if( cy < y_min ) { s = cy - y_min; d += s * s; } else
      if( cy > y_max ) { s = cy - y_max; d += s * s; }

      if( cz < z_min ) { s = cz - z_min; d += s * s; } else
      if( cz > z_max ) { s = cz - z_max; d += s * s; }

      return d <= r * r;
    }
    bool sphere_intersect(float r, float cx, float cy, float cz) {
      AABB bbox(cx - r, cx + r, 
		cy - r, cy + r,
		cz - r, cz + r);
      return sphere_intersect(bbox, r, cx, cy, cz);
    }


    // Intersect two AABB. 
    bool intersect_AABB(AABB &R) {
      // Negation of (not intersection). 
      return !( R.x_max < x_min || R.x_min > x_max ) && 
	     !( R.y_max < y_min || R.y_min > y_max ) &&
	     !( R.z_max < z_min || R.z_min > z_max ); 
    }

    // Returns true if the current region encloses the passed region R. 
    bool encloses_AABB(AABB &R ){ 
      return 
	( (R.x_max <= x_max) && (R.x_min >= x_min) ) && 
	( (R.y_max <= y_max) && (R.y_min >= y_min) ) &&
	( (R.z_max <= z_max) && (R.z_min >= z_min) );
    }

    // Assumes caller computed a minimum AABB around the triangle as
    // follows: AABB tribox(tri.min_x(), tri.max_x(), tri.min_y(),
    // tri.max_y(), tri.min_z(), tri.max_z());
    bool tri_intersect(triangle &tri, AABB &tribox) {
      // 1. Does the AABB enclose the minimum AABB triangle box. 
      if(encloses_AABB(tribox)) return true;

      // 2. Test overlap of minimum AABB of the triangle.
      if(intersect_AABB(tribox) == false) return false; // No overlap, no intersection.

      // 3. Apply full Moller intersection testing
      float boxcenter[3], boxhalfsize[3], triverts[3][3];
      boxhalfsize[0] = (x_max - x_min) / 2.0;
      boxhalfsize[1] = (y_max - y_min) / 2.0;
      boxhalfsize[2] = (z_max - z_min) / 2.0;
      boxcenter[0] = x_min + boxhalfsize[0];
      boxcenter[1] = y_min + boxhalfsize[1];
      boxcenter[2] = z_min + boxhalfsize[2];
      for(int i = 0; i < 3; i++) 
	for(int j = 0; j < 3; j++) triverts[i][j] = tri[i][j];


      // TODO: We can get rid of triBoxOverlap and tribox3.cpp using
      // this alternative algorithm.  This is *SIMPLER* probably 2x
      // slower.  Create a small triangle mesh of the AABB. This will
      // have 2 x 6 = 12 triangles.  Cast a ray from the 3 edges of
      // the intersecting triangle.  If the edge ray intersects the
      // triangle with distance < the edge length. Then there is an
      // intersection. We can also ask if the minimum box of those
      // triangles intersects the min box of the intersecting triangle

      // TODO: Can also replace with a 6 x 2 triangle/triangle
      // intersections.
      // TODO: Find a simple triangle/triangle intersection algorithm.

      return triBoxOverlap(boxcenter, boxhalfsize, triverts) > 0;
    }       
  };

  // T must support x(), y(), z(), min/max_x(), min/max_y(), min/max_z().
  template <typename T> struct BVH {
    // Comparitors for partitioning data x,y, or z axes.
    struct x_cmp { bool operator () (T a, T b) const { return a.x() < b.x(); } };
    struct y_cmp { bool operator () (T a, T b) const { return a.y() < b.y(); } };
    struct z_cmp { bool operator () (T a, T b) const { return a.z() < b.z(); } };
    typedef typename vector<T>::iterator it_t;
    
    vector<T> &data;

    // Lmis, Rmids, Lbv, and Rbv all implement a linear probing hash
    // table.  Given a midpoint the hash table returns an axis aligned
    // bounding box.
    vector<it_t> Lmids, Rmids; vector< AABB > Lbv, Rbv;
    size_t leftN, rightN;

    void insert(vector<it_t> &mids, vector<AABB> &bv,
		it_t mid, AABB &aabb) {
      // Get the distance of th enew midpoint from the start of the data array.
      size_t i = distance(data.begin(), mid) % mids.size();
      // Find an open slot.
      size_t istart = i;
      while(mids[i] != data.end()) { 
	i = (i+1) % mids.size();
	if(i == istart) { cerr << "BVH: split table is full!" << endl; exit(1); }
      }
      // Save midpoint and bounding box.
      mids[i] = mid; bv[i] = aabb;
    }

    // Look for AABB for a given mid point in hash table.
    AABB find(vector<it_t> &mids, vector<AABB> &bv, it_t mid) {
      size_t i = distance(data.begin(), mid) % mids.size();
      size_t istart = i;
      while(mids[i] != data.end()) {
	if(mids[i] == mid) return bv[i]; else i = (i + 1) % mids.size();
	if(i == istart) { cerr << "BVH: Did not find mid. Problem!!!" << endl; exit(1);}
      }
      return AABB();
    }

    void build(it_t first, it_t last, int dim) {
      if (last - first <= 7) return; // Since 3 dimensions are partitioned.

      it_t mid = first + ( last - first + 1 ) / 2;
      // NOTE: nth_elment calls the destructor if it exists. This
      // sometimes creates issues.
      switch(dim) { // Cycle through x, y, z axes.
      case 0: nth_element(first, mid, last, x_cmp()); dim++;   break;
      case 1: nth_element(first, mid, last, y_cmp()); dim++;   break;
      case 2: nth_element(first, mid, last, z_cmp()); dim = 0; break;
      }

      // For a given split point mid, save in a hash table the right and left AABBs.
      AABB left;  BV(left, first, mid); insert(Lmids, Lbv, mid, left);  leftN++;
      AABB right; BV(right, mid, last); insert(Rmids, Rbv, mid, right); rightN++;

      build(first, mid, dim); build(mid, last, dim);
    }

    // Determine bounding volume, AABB given in iterator range.
    void BV(AABB &bv, it_t first, it_t last) {
      // Maximum bounding range.
      bv.x_max = -(numeric_limits<float>::max()); bv.x_min = numeric_limits<float>::max();
      bv.y_max = -(numeric_limits<float>::max()); bv.y_min = numeric_limits<float>::max();
      bv.z_max = -(numeric_limits<float>::max()); bv.z_min = numeric_limits<float>::max();

      // Scan data to determine dimensions of AABB.
      for(it_t it = first; it != last; it++) {
	T &obj = *it;
	float min_x = obj.min_x(), min_y = obj.min_y(), min_z = obj.min_z();
	if(min_x < bv.x_min) bv.x_min = min_x; 
	if(min_y < bv.y_min) bv.y_min = min_y; 
	if(min_z < bv.z_min) bv.z_min = min_z;

	float max_x = obj.max_x(), max_y = obj.max_y(), max_z = obj.max_z();
	if(max_x > bv.x_max) bv.x_max = max_x; 
	if(max_y > bv.y_max) bv.y_max = max_y; 
	if(max_z > bv.z_max) bv.z_max = max_z;
      }
    }

    BVH(vector<T> &data_) : data(data_) {
      leftN = rightN = 0;
      
      Lbv.resize(data.size()/2); Lmids.resize(data.size()/2, data.end());
      Rbv.resize(data.size()/2); Rmids.resize(data.size()/2, data.end());

      build(data.begin(), data.end(), 0);
    }

    // Recursively intersect a ray with members of the BVH.
    void query(it_t first, it_t last,  
	       float x, float y, float z, // origin
	       float idx, float idy, float idz, // direction with each dimensio inverted.
	       vector<T*> &res)  {
      if( last - first <= 7 ) { 
	for(it_t it = first; it != last; it++) {
	  AABB test(it->min_x(), it->max_x(), it->min_y(), it->max_y(), it->min_z(), it->max_z());
	  if(test.ray_intersect(x,y,z,idx,idy,idz)) res.push_back(&*it);  
	}
	return; 
      }

      it_t mid = first + ( last - first + 1 ) / 2;
      // Given a mid point "mid" find an axis aligned bounding box for the "left" side.
      AABB left = find(Lmids, Lbv, mid);
      if(left.ray_intersect(x,y,z,idx,idy,idz))  query(first, mid, x,y,z, idx,idy,idz, res);

      // Given a mid point "mid" find an axis aligned bounding box for the "right" side.
      AABB right = find(Rmids, Rbv, mid);
      if(right.ray_intersect(x,y,z,idx,idy,idz)) query(mid, last, x,y,z, idx,idy,idz, res);
    }

    // Recursively intersect a sphere with BVH members.
    void query(it_t first, it_t last,  
	       float r, 
	       float cx, float cy, float cz,
	       vector<T*> &res)  {
      if( last - first <= 7 ) { 
	for(it_t it = first; it != last; it++) {
	  T &v = *it;
	  AABB vbox(v.min_x(), v.max_x(), v.min_y(), v.max_y(), v.min_z(), v.max_z());
	  if(vbox.sphere_intersect(r,cx,cy,cz))  res.push_back(&v);  
	}
	return; 
      }
      it_t mid = first + ( last - first + 1 ) / 2;
      AABB left = find(Lmids, Lbv, mid);
      if(left.sphere_intersect(r,cx,cy,cz))  query(first, mid, r, cx,cy,cz , res);
      AABB right = find(Rmids, Rbv, mid);
      if(right.sphere_intersect(r,cx,cy,cz)) query(mid, last, r, cx,cy,cz , res);
    }

    // Cast a ray through bounding volume with origin [x,y,z]' and direction [dx,dy,dz]'.
    void Query(vector<T*> &res,
	       float x, float y, float z,    // origin 
	       float dx, float dy, float dz) { // direction
      query(data.begin(), data.end(), x,y,z, 1.0f/dx, 1.0f/dy, 1.0f/dz, res); 
    }

    // Intersects sphere with AABBs of the BVH.
    void Query(vector<T*> &res, float r,  // radius of sphere
	       float cx, float cy, float cz) { // center of sphere
      query(data.begin(), data.end(), r, cx,cy,cz, res); 
    }
  };

  // Pointer based interface for a BVH. 
  template <typename T> struct BVHptr {
    T *ptr;
    // Accessors for BVH<> 
    float x() { return ptr->x(); } 
    float y() { return ptr->y(); } 
    float z() { return ptr->z(); }
    float min_x() { return ptr->min_x(); } float max_x() { return ptr->max_x(); }
    float min_y() { return ptr->min_y(); } float max_y() { return ptr->max_y(); }
    float min_z() { return ptr->min_z(); } float max_z() { return ptr->max_z(); }        
  };

  // BVH2 allows pointers to be used for a BVH. This is accomplished
  // by "wrapping" th epointer in the BVHptr class above.
  template <typename T> struct BVH2 { 
    vector< BVHptr<T> > ptrs;
    BVH< BVHptr<T> > *bvh;
    
    BVH2(vector <T> &data) {
      ptrs.resize(data.size());
      for(size_t i = 0; i < ptrs.size(); i++) ptrs[i].ptr = &data[i];
      bvh = new BVH< BVHptr<T> > (ptrs);
    }
    BVH2(vector <T *> &data) {
      ptrs.resize(data.size());
      for(size_t i = 0; i < ptrs.size(); i++) ptrs[i].ptr = data[i];
      bvh = new BVH< BVHptr<T> > (ptrs);
    }
    ~BVH2() { delete bvh; }

    void Query(vector<T*> &res,  
	       float x, float y, float z,    // origin 
	       float dx, float dy, float dz) { // direction
      vector< BVHptr<T> * > res_q;
      bvh->Query(res_q, x, y, z, dx, dy, dz);
      for(size_t i = 0; i < res_q.size(); i++) res.push_back(res_q[i]->ptr);
    }

    void Query(vector<T*> &res, float r,  // radius of sphere
	       float cx, float cy, float cz) { // center of sphere
      vector< BVHptr<T> * > res_q;
      bvh->Query(res_q, r, cx, cy, cz);
      for(size_t i = 0; i < res_q.size(); i++) res.push_back(res_q[i]->ptr);
    }
	       
  };
}

#endif
