#include <vector>
#include <iomanip>

#include "geometry.hpp"

namespace geom {
  using namespace std;

  // 3x3 matrix inverse
  matrix3 inv(const matrix3 &A) {
    float detA = A.det();
    float a11 = A(0,0), a12 = A(0,1), a13 = A(0,2);
    float a21 = A(1,0), a22 = A(1,1), a23 = A(1,2);
    float a31 = A(2,0), a32 = A(2,1), a33 = A(2,2);

    matrix3 R;
    R(0,0) = det2(a22, a23, a32, a33);  R(0,1) = det2(a13, a12, a33, a32);  R(0,2) = det2(a12, a13, a22, a23);
    R(1,0) = det2(a23, a21, a33, a31);  R(1,1) = det2(a11, a13, a31, a33);  R(1,2) = det2(a13, a11, a23, a21);  
    R(2,0) = det2(a21, a22, a31, a32);  R(2,1) = det2(a12, a11, a32, a31);  R(2,2) = det2(a11, a12, a21, a22);

    return (1.0/detA) * R;
  }

  // Compute coefficients of characteristic polynomial for a 3x3 matrix. NOTE: a = -1
  void CharPoly3(matrix3 &A, float &b, float &c, float &d) {
    float trA = A.trace();
    matrix3 AxA; AxA = A * A;
    b = trA; c = 0.5 * ( AxA.trace() - trA * trA);  d = A.det();
  }

  // Returns false if roots are complex. They shouldn't be for a
  // covariance matrix.  Roots are sorted from smallest to largest.
  bool SolveCubic(vector<float> &x, float a, float b, float c, float d) {
    float f = (3 * c / a - b*b / (a*a)) / 3;
    float g = (2 * b*b*b  / (a*a*a) - 9 * b * c / (a*a) + 27 * d / a) / 27;
    float h = g*g / 4 + f*f*f / 27;
    if(h <= 0) {
      float i = sqrt(g*g / 4 - h);
      float j = pow(i, 1.0 / 3.0);
      float K = acos(-g / (2*i));
      float L = j * - 1;
      float M = cos(K/3);
      float N = sqrt(3) * sin(K/3);
      float P = b / (3*a) * -1;
      x.resize(3);
      x[0] = 2*j*cos(K/3) - b/(3*a); x[1] = L * (M+N) + P; x[2] = L * (M-N) + P; // Roots..
      sort(x.begin(),x.end());
      return true;
    }
    else return false;
  }

  // Given a guess for an eigenvalue mu, compute the corresponding eigenvector.
  vec3 InverseIteration(matrix3 &A, float mu) {
    vec3 b(0.3245,0.5432,0.1234);
    matrix3 Q = inv(A - mu * (eye3()));
    b = Q * b; b.normalize(); // iteration 1
    b = Q * b; b.normalize(); // iteration 2
    b = Q * b; b.normalize(); // iteration 2
    return b;
  }

  // Compute the covariance of the given set of 3-d vectors. 
  matrix3 cov(vector<vec3> &data) {
    vec3 mu;
    float N = data.size();
    // Compute mean.
    for(size_t i = 0; i < data.size(); i++) mu += data[i] / N;
    matrix3 R; // Compute covariance.
    for(int i = 0; i < 3; i++) {
      for(int j = i; j < 3; j++) {
	for(size_t k = 0; k < data.size(); k++) {
	  R(i,j) += ((data[k][i] - mu[i]) * (data[k][j] - mu[j])) / N;
	}
	R(j,i) = R(i,j); // Symmetric matrix.
      }
    }
    return R;
  }

  vec3 Centroid3(vector<vec3> &points) {
    vec3 C(0,0,0); float N = points.size();
    for(size_t i = 0; i < points.size(); i++) C += points[i] / N;
    return C;
  }

  vec3 Centroid3(vector<ivec3> &points) {
    vec3 C(0,0,0); float N = points.size();
    for(size_t i = 0; i < points.size(); i++) C += points[i].to_vec3() / N;
    return C;
  } 

  void ClearLargest(vector< vector<ivec3> > &comps) {
    if(comps.size() == 0) return;
    size_t max_sz = comps[0].size();
    for(size_t c = 1; c < comps.size(); c++) max_sz = max(max_sz, comps[c].size());
    for(size_t c = 0; c < comps.size(); c++) {
      if(comps[c].size() == max_sz) comps[c].clear();
    }
  }

  bool PCA3(vector<vec3> &v, vector<float> &x, vector<vec3> &points) {
    matrix3 C = cov(points);
    v.resize(3); x.resize(3);
    // Compute characteristic polyonomial of 3x3 covariance matrix. 
    float b,c,d;  CharPoly3(C, b,c,d); // NOTE: a = -1
    if(SolveCubic(x, -1, b,c,d)) {
      // Perturb eigenvalues.
      vector<float> mu(3);
      const float alpha = 0.00001;
      mu[0] = (1 - alpha) * x[0] + alpha * x[1];
      mu[1] = (1 - alpha) * x[1] + alpha * x[2];
      mu[2] = alpha * x[1] + (1 - alpha) * x[2];
      // Get corresponding eigenvector by inverse iteration.
      for(int i = 0; i < 3 ; i++) v[i] = InverseIteration(C, mu[i]); 
      return true;
    }
    else return false;
  }

  matrix3 outer3(vec3 &a, vec3 &b) {
    matrix3 A;
    A(0,0) = a[0] * b[0]; A(0,1) = a[0] * b[1]; A(0,2) = a[0] * b[2];
    A(1,0) = a[1] * b[0]; A(1,1) = a[1] * b[1]; A(1,2) = a[1] * b[2];
    A(2,0) = a[2] * b[0]; A(2,1) = a[2] * b[1]; A(2,2) = a[2] * b[2];
    return A;
  }

  // Area of a triangle in 3d.
  float area(vec3 &a, vec3 &b, vec3 &c) {
    float xA = a.x(), yA = a.y(), zA = a.z();
    float xB = b.x(), yB = b.y(), zB = b.z();
    float xC = c.x(), yC = c.y(), zC = c.z();
    float det0 = det3(xA, xB, xC,
		       yA, yB, yC,
		       1 ,  1,  1);
    float det1 = det3(yA, yB, yC,
		       zA, zB, zC,
		       1 ,  1,  1);
    float det2 = det3(zA, zB, zC, 
		       xA, xB, xC,
		       1 ,  1,  1);
    return 0.5f * sqrtf( det0 * det0 + det1 * det1 + det2 * det2 );
  }

};

// Implements an in-pace 3d-tree for 3-d point data. Supports fast
// orthogonal range queries.
namespace spacepartition {
  /*
    < class T > must implement the methods
    float x() const { return xval; }  float y() const { return yval; }  float z() const { return zval; }
  */

  using namespace std;

  enum dim {x_axis, y_axis, z_axis}; // Axis on which to conduct the split.
  // Specifies the side of the spliting region left is less than right is greater than.
  enum region_t { right_region, left_region };

  // 3-d orthogonal region
  template <class T> struct region {
    float x_min, x_max, y_min, y_max, z_min, z_max;

    // Simplest constructor. It just assigns max and min elements as specified.
    region(float x_min_, float x_max_, 
	   float y_min_, float y_max_,
	   float z_min_, float z_max_){ 
      x_min = x_min_; x_max = x_max_; y_min = y_min_; y_max = y_max_; z_min = z_min_; z_max = z_max_;
    }
    // If no ranges are specified, the region spans the entire (x, y) space.
    region() { 
      x_min = -(numeric_limits<float>::max()); x_max = numeric_limits<float>::max();
      y_min = -(numeric_limits<float>::max()); y_max = numeric_limits<float>::max();
      z_min = -(numeric_limits<float>::max()); z_max = numeric_limits<float>::max();
    }

    // Create a half plane region left, right, top or bottom.
    region(float splitline, dim splitdim, region_t r) {
      x_min = -(numeric_limits<float>::max()); x_max = numeric_limits<float>::max();
      y_min = -(numeric_limits<float>::max()); y_max = numeric_limits<float>::max();
      z_min = -(numeric_limits<float>::max()); z_max = numeric_limits<float>::max();

      // Figure out which axis the node splits on.
      if(splitdim == x_axis){
	// Right if a "right" region region (greater than or equal to).
	if(r == right_region) x_min = splitline;
	else if( r == left_region )  x_max = splitline;
	// Left if a "left" region (strictly less than). 
      }
      else if(splitdim == y_axis){
	// Top if a "right" region region (greater than or equal to).
	if(r == right_region) y_min = splitline;
	else if( r == left_region )  y_max = splitline;
	// Bottom if a "left" region (strictly less than). 
      }
      else if(splitdim == z_axis) {
	if(r == right_region) z_min = splitline;
	else if( r == left_region )  z_max = splitline;
      }
    }

    // Returns true if the current region encloses the passed region R. 
    bool encloses( region R ){ 
      return 
	( (R.x_max <= x_max) && (R.x_min >= x_min) ) && 
	( (R.y_max <= y_max) && (R.y_min >= y_min) ) &&
	( (R.z_max <= z_max) && (R.z_min >= z_min) );
    }

    // Returns true if the current region contains a point (or pointer to a point).  
    bool encloses( T p ) { return 
	( (p.x()  <= x_max) && (p.x()  >= x_min) ) && 
	( (p.y()  <= y_max) && (p.y()  >= y_min) ) &&
	( (p.z()  <= z_max) && (p.z()  >= z_min) ); 
    }
    bool encloses( T *p ){ return 
	( (p->x() <= x_max) && (p->x() >= x_min) ) && 
	( (p->y() <= y_max) && (p->y() >= y_min) ) &&
	( (p->z() <= z_max) && (p->z() >= z_min) ); 
    }

    // Returns true of the current region and R intersect. 
    bool intersects( region &R ){    
      // Negation of (not intersection). 
      return 
	!( R.x_max < x_min || R.x_min > x_max) && 
	!( R.y_max < y_min || R.y_min > y_max ) &&
	!( R.z_max < z_min || R.z_min > z_max ); 
    }
    // Constructor creates a new region that is the intersection of the two given regions.
    region ( region a, region b) {
      x_max = min(a.x_max, b.x_max); x_min = max(a.x_min, b.x_min);
      y_max = min(a.y_max, b.y_max); y_min = max(a.y_min, b.y_min);
      z_max = min(a.z_max, b.z_max); z_min = max(a.z_min, b.z_min);
    }
  };

  template <class T> struct _3dtree {
    struct x_cmp { bool operator () (T a, T b) const { return a.x() < b.x(); } };
    struct y_cmp { bool operator () (T a, T b) const { return a.y() < b.y(); } };
    struct z_cmp { bool operator () (T a, T b) const { return a.z() < b.z(); } };
    typedef typename vector<T>::iterator it_t;

    void build(it_t first, it_t last) {
      if(last - first <= 7) return; // <= 7, since need mid, Lmid, and Rmid

      it_t mid  = first + ( last - first + 1 ) / 2;
      it_t Lmid = first + ( mid - first + 1 ) / 2;  // first to mid
      it_t LLmid = first + ( Lmid - first + 1) / 2; // first to Lmid
      it_t RLmid = Lmid + (mid - Lmid + 1) / 2;     // Lmid to mid
      it_t Rmid = mid + ( last - mid + 1 ) / 2;     // mid to last
      it_t LRmid = mid + ( Rmid - mid + 1) / 2;;    // mid to Rmid
      it_t RRmid = Rmid + ( last - Rmid + 1) / 2;   // Rmid to last
    
      nth_element(first, mid, last,  x_cmp());

      nth_element(first, Lmid, mid,  y_cmp());   // first to mid
      nth_element(first, LLmid, Lmid, z_cmp());  // first to Lmid
      nth_element(Lmid+1, RLmid, mid, z_cmp());  // Lmid to mid

      nth_element(mid+1, Rmid, last, y_cmp());   // mid to last
      nth_element(mid+1, LRmid, Rmid, z_cmp());  // mid to Rmid
      nth_element(Rmid+1, RRmid, last, z_cmp()); // Rmid to last

      build(first, LLmid); build(LLmid+1, Lmid); build(Lmid+1, RLmid); build(RLmid+1, mid);
      build(mid+1, LRmid); build(LRmid+1, Rmid); build(Rmid+1, RRmid); build(RRmid+1, last);
    }

    vector<T> &P;
    region<T> region0;

    _3dtree(vector<T> &P_) : P(P_) {   
      region0.x_max = -(numeric_limits<float>::max()); region0.x_min = numeric_limits<float>::max();
      region0.y_max = -(numeric_limits<float>::max()); region0.y_min = numeric_limits<float>::max();
      region0.z_max = -(numeric_limits<float>::max()); region0.z_min = numeric_limits<float>::max();

      for(size_t i = 0; i < P.size(); i++) {
	float x = P[i].x(), y = P[i].y(), z = P[i].z();
	if(x < region0.x_min) region0.x_min = x; if(y < region0.y_min) region0.y_min = y; if(z < region0.z_min) region0.z_min = z;
	if(x > region0.x_max) region0.x_max = x; if(y > region0.y_max) region0.y_max = y; if(z > region0.z_max) region0.z_max = z;
      }
      // cout << "_3dtree(build): P.size()=" << P.size() << endl;
      build(P.begin(), P.end());   
    }

    // Orothogonal range query.
    void query(it_t first, it_t last, region<T> vr, region<T> R, vector<T*> &res)  {
      if(R.encloses(vr)) { // Region is entirely enclosed, report everyone.
	for(it_t it = first; it != last; it++) res.push_back(&*it);
	return;
      }
      if( last - first <= 7 ) {
	for(it_t it = first; it != last; it++) if(R.encloses(*it)) res.push_back(&*it);
	return;
      }

      it_t mid  = first + ( last - first + 1 ) / 2;
      it_t Lmid = first + ( mid - first + 1 ) / 2;  // first to mid
      it_t LLmid = first + ( Lmid - first + 1) / 2; // first to Lmid
      it_t RLmid = Lmid + (mid - Lmid + 1) / 2;     // Lmid to mid
      it_t Rmid = mid + ( last - mid + 1 ) / 2;     // mid to last
      it_t LRmid = mid + ( Rmid - mid + 1) / 2;;    // mid to Rmid
      it_t RRmid = Rmid + ( last - Rmid + 1) / 2;   // Rmid to last

      // Create regions based on specified splitting planes.
      region<T> left( region<T>((*mid).x(), x_axis, left_region),  vr);
      region<T> right(region<T>((*mid).x(), x_axis, right_region), vr);
      region<T> Lbot(region<T>((*Lmid).y(), y_axis, left_region), left);      // 1
      region<T> Lbotbot(region<T>((*LLmid).z(), z_axis, left_region), Lbot);  // LLmid
      region<T> Lbottop(region<T>((*LLmid).z(), z_axis, right_region), Lbot); // 2
      region<T> Ltop(region<T>((*Lmid).y(), y_axis, right_region), left);     // 3
      region<T> Ltopbot(region<T>((*RLmid).z(), z_axis, left_region), Ltop);  // RLmid
      region<T> Ltoptop(region<T>((*RLmid).z(), z_axis, right_region), Ltop); // 4

      region<T> Rbot(region<T>((*Rmid).y(), y_axis, left_region), right);     // 5
      region<T> Rbotbot(region<T>((*LRmid).z(), z_axis, left_region), Rbot);  // LRmid
      region<T> Rbottop(region<T>((*LRmid).z(), z_axis, right_region), Rbot); // 6
      region<T> Rtop(region<T>((*Rmid).y(), y_axis, right_region), right);    // 7
      region<T> Rtopbot(region<T>((*RRmid).z(), z_axis, left_region), Rtop);  // RRmid
      region<T> Rtoptop(region<T>((*RRmid).z(), z_axis, right_region), Rtop); // 8

      // 8-way recursion for 3-d data.
      if(R.intersects(Lbotbot)) query(first, LLmid, Lbotbot, R, res);
      if(R.intersects(Lbottop)) query(LLmid+1, Lmid, Lbottop, R, res);
      if(R.intersects(Ltopbot)) query(Lmid+1, RLmid, Ltopbot, R, res);
      if(R.intersects(Ltoptop)) query(RLmid+1, mid, Ltoptop, R, res);

      if(R.intersects(Rbotbot)) query(mid+1, LRmid, Rbotbot, R, res);
      if(R.intersects(Rbottop)) query(LRmid+1, Rmid, Rbottop, R, res);
      if(R.intersects(Rtopbot)) query(Rmid+1, RRmid, Rtopbot, R, res);
      if(R.intersects(Rtoptop)) query(RRmid+1, last, Rtoptop, R, res);

      if(R.encloses(*mid))  res.push_back(&*mid);
      if(R.encloses(*Lmid)) res.push_back(&*Lmid);
      if(R.encloses(*LLmid)) res.push_back(&*LLmid);
      if(R.encloses(*RLmid)) res.push_back(&*RLmid);
      if(R.encloses(*Rmid)) res.push_back(&*Rmid);
      if(R.encloses(*LRmid)) res.push_back(&*LRmid);
      if(R.encloses(*RRmid)) res.push_back(&*RRmid);
    }

    // Performs an orthogonal range query.
    void Query ( region< T > R, vector<T*> &res) {  query(P.begin(), P.end(), region0, R, res);   }
    // Brute force orthogonal range query.
    void SlowQuery( region< T > R, vector<T*> &res) {
      for(size_t i = 0; i < P.size(); i++) if(R.encloses(P[i])) res.push_back(&P[i]);
    }
  };


  //void F(vector<float> &/*dydt*/, float /* t */, vector<float> &/*y*/) { }

  // Adapted from ODE Solver: RK45 version 1.1.  Connelly Barnes 2005.
  // Downloaded from http://www.connellybarnes.com/code/python/rk45
  // Sep 22nd, 2010.  Also adapted methods from
  // http://www.library.cornell.edu/nr/bookcpdf/c16-2.pdf
  /*void rk45(vector<float> &y, vector<float> &y0, float t0, float tfinal, float tol) {
    y = y0;
    int i = 0, N = y0.size();
    float t = t0;
    float hmax = (tfinal - t0) / 128.0;
    float h = hmax / 4.0;
    // Cash-Karp parameters
    const float a[]  = { 0.0, 0.2, 0.3, 0.6, 1.0, 0.875 };
    const float b1[] = {0.2};
    const float b2[] = {3.0/40.0, 9.0/40.0};
    const float b3[] = {0.3, -0.9, 1.2};
    const float b4[] = {-11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0};
    const float b5[] = {1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0};
    const float c []  = {37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0};
    const float dc[] = {c[0]-2825.0/27648.0, c[1]-0.0, c[2]-18575.0/48384.0,
 			c[3]-13525.0/55296.0, c[4]-277.00/14336.0, c[5]-0.25};
 
    vector< vector<float> > k(5); for(i = 0; i < 5; i++) k[i].resize(N);
    vector<float> tmp(N);
 
    while(t < tfinal) {
      if (t + h > tfinal) h = tfinal - t;
      if (t + h <= t) { cerr << "Singularity in ODE" << endl; exit(1); }
       
      // k[0] = F(t, y)
      F(k[0], t, y); 
       
      // k[1] = F(t+a[1]*h, y+h*(k[0]*b[1][0]))
      for(i = 0; i < N; i++) tmp[i] = y[i] + h * k[0][i] * b1[0]; 
      F(k[1], t+a[1]*h, tmp);
 
      // k[2] = F(t+a[2]*h, y+h*(k[0]*b[2][0]+k[1]*b[2][1]))
      for(i = 0; i < N; i++) tmp[i] = y[i]+h*(k[0][i]*b2[0]+k[1][i]*b2[1]);
      F(k[2], t+a[2]*h, tmp);
 
      // k[3] = F(t+a[3]*h, y+h*(k[0]*b[3][0]+k[1]*b[3][1]+k[2]*b[3][2]))
      for(i = 0; i < N; i++) tmp[i] = y[i]+h*(k[0][i]*b3[0]+k[1][i]*b3[1]+k[2][i]*b3[2]);
      F(k[3], t+a[3]*h, tmp);
 
      // k[4] = F(t+a[4]*h, y+h*(k[0]*b[4][0]+k[1]*b[4][1]+k[2]*b[4][2]+k[3]*b[4][3]))      
      for(i = 0; i < N; i++) tmp[i] = y[i]+h*(k[0][i]*b4[0]+k[1][i]*b4[1]+k[2][i]*b4[2]+k[3][i]*b4[3]);
      F(k[4], t+a[4]*h, tmp);
 
      // k[5] = F(t+a[5]*h, y+h*(k[0]*b[5][0]+k[1]*b[5][1]+k[2]*b[5][2]+k[3]*b[5][3]+k[4]*b[5][4]))
      for(i = 0; i < N; i++) tmp[i] = y[i]+h*(k[0][i]*b5[0]+k[1][i]*b5[1]+k[2][i]*b5[2]+k[3][i]*b5[3]+k[4][i]*b5[4]);
      F(k[5], t+a[5]*h, tmp);
 
      // Estimate current error and current maximum error.
      // E = norm(h*(k[0]*dc[0]+k[1]*dc[1]+k[2]*dc[2]+k[3]*dc[3]+k[4]*dc[4]+k[5]*dc[5]))      
      for(i = 0; i < N; i++) tmp[i] = h*(k[0][i]*dc[0]+k[1][i]*dc[1]+k[2][i]*dc[2]+k[3][i]*dc[3]+k[4][i]*dc[4]+k[5][i]*dc[5]);
      float E = 0; for(i = 0; i < N; i++) E += tmp[i] * tmp[i]; E = sqrt(E);
 
      // Emax = tol*max(norm(y), 1.0)      
      float normy = 0; for( i = 0; i < N; i++) normy += y[i] * y[i]; normy = sqrt(normy);
      float Emax = tol*max(normy, (float)1.0);
 
      // Update solution if error is OK.
      if (E < Emax) {
 	t += h;
 	// y += h*(k[0]*c[0]+k[1]*c[1]+k[2]*c[2]+k[3]*c[3]+k[4]*c[4]+k[5]*c[5])
 	for(i = 0; i < N; i++) y[i] += h*(k[0][i]*c[0]+k[1][i]*c[1]+k[2][i]*c[2]+k[3][i]*c[3]+k[4][i]*c[4]+k[5][i]*c[5]);
      } 
      // Update step size
      if (E > 0.0) { h = min(hmax, float(0.85*h*pow((Emax/E),0.2))); }
    }
    }*/

}
