#include "procthreads.hpp"

#include <QThread>
// NOTE: On Mac OSX, cout/cin are not thread safe. They will cause
// crashes if you output in multithreaded code.

// Given a subset of cells an a refinement volume16, applied image
// refinement to cell.
struct RefineThread : public QThread {
  vector<mesh *> cells;
  void addCell(mesh *cell) { cells.push_back(cell); }
  volume8 &v;
  float dist; int max_iter; float alpha;
  float stepsize;
  RefineThread(volume8 &v_, float dist, int max_iter, float alpha, float stepsize) : 
    v(v_), dist(dist), max_iter(max_iter), alpha(alpha), stepsize(stepsize)  {  }
  void run() { for(size_t i = 0; i < cells.size(); i++) cells[i]->refine(v, dist, max_iter, alpha, stepsize); }
};

void refine(volume8 *v, float dist, int max_iter, float alpha, float stepsize, vector<geom::mesh *> &cells, int num_threads) {
  if(cells.size() == 0) return;
  vector<RefineThread *> threads(num_threads);
  for(int T = 0; T < num_threads; T++) threads[T] = new RefineThread(*v, dist, max_iter, alpha, stepsize);
  int T = 0;
  for(size_t i = 0; i < cells.size(); i++) {
    threads[T]->addCell(cells[i]);
    T++; if(T >= num_threads) T = 0;
  }
  for(int t = 0; t < num_threads; t++) threads[t]->start();
  for(int t = 0; t < num_threads; t++) threads[t]->wait();
  for(int t = 0; t < num_threads; t++) delete threads[t]; 
};

// Applies histogram contrast enhancement to each step scan line.
struct HistApplyThread : public QThread {
  volume8 *v; 
  histvol_t<uint8> &hv;
  int startScanLine, stepScanLine;  
  HistApplyThread(int start, int step, volume8 *v_, histvol_t<uint8> &hv_) : hv(hv_) { 
    v = v_; startScanLine = start; stepScanLine = step;
  }
  void run() { v->ApplyHistogram(hv, startScanLine, stepScanLine); }
};

void happly(volume8 *v, int nx, int ny, int nz, float fmin, float fmax, int num_threads) {
  // Make sure intensity histogram is populated first!
  histvol_t<uint8> hv(nx, ny, nz, fmin, fmax);
  v->PopulateHistVol(hv);
  hv.pre_process(); 
  vector<HistApplyThread *> threads(num_threads);
  for(int t = 0; t < num_threads; t++) threads[t] = new HistApplyThread(t, num_threads, v, hv);
  for(int t = 0; t < num_threads; t++) threads[t]->start();
  for(int t = 0; t < num_threads; t++) threads[t]->wait();
  for(int t = 0; t < num_threads; t++) delete threads[t];
}


struct StretchThread : public QThread {
  histvol_t<uint16> &hv;
  volume16 *v; 
  int startLine, stepLine;  
  StretchThread(int startLine, int stepLine, volume16 *v, histvol_t<uint16> &hv) :  
    hv(hv), v(v), startLine(startLine), stepLine(stepLine) { }
  void run() { v->ApplyHistogram(hv, startLine, stepLine); }  
};

void stretch16(volume16 *v, int num_threads) {
  cout << "*** stretch16(): performing *NOW* ***" << endl;
  histvol_t<uint16> hv(1, 1, 1, 0.001, 0.001);
  v->PopulateHistVol(hv);
  hv.pre_process();
  cout << "hv.vmin=" << hv(0,0,0).vmin << " hv.vmax=" << hv(0,0,0).vmax << endl;
  vector<StretchThread *> threads(num_threads);
  for(int t = 0; t < num_threads; t++) threads[t] = new StretchThread(t, num_threads, v, hv);
  for(int t = 0; t < num_threads; t++) threads[t]->start();
  for(int t = 0; t < num_threads; t++) threads[t]->wait();
  for(int t = 0; t < num_threads; t++) delete threads[t];
}


// Thread applies Gaussian filter to every step scan line.
struct GaussThread : public QThread {
  volume8 *v;
  gauss_dim dim;
  int startScanLine, stepScanLine;  
  float sigma;
  GaussThread(gauss_dim dim_, int start, int step, volume8 *v_, float sigma_) {
    v = v_; sigma = sigma_; dim = dim_; startScanLine = start; stepScanLine = step;
  }
  void run() {
    switch(dim) {
    case DimX: v->GaussX(sigma, startScanLine, stepScanLine); break;
    case DimY: v->GaussY(sigma, startScanLine, stepScanLine); break;
    case DimZ: v->GaussZ(sigma, startScanLine, stepScanLine); break;
    default: break;
    }
  }
};

void blur(volume8 *v, gauss_dim dim, float sigma, int num_threads) {
  vector<GaussThread *> threads;
  for(int tid = 0; tid < num_threads; tid++) threads.push_back(new GaussThread(dim, tid, num_threads, v, sigma));
  cout << "blur threads=" << threads.size() << endl;
  for(int t = 0; t < num_threads; t++) threads[t]->start();
  for(int t = 0; t < num_threads; t++) threads[t]->wait();
  for(int t = 0; t < num_threads; t++) delete threads[t];
}


// Performs dilate on every step scan line.
class DilateThread : public QThread {
  int zstart, zstep, r, bg, fg; float alpha; volume8 &v; bool invert;
  void run() { v.dilate(res, r, fg, bg, alpha, invert, zstart, zstep); }
public:
  vector<ivec3> res;
  DilateThread(volume8 &vi,
	       int r, int bg, int fg, float alpha, bool invert, int zstart, int zstep) : v(vi) {
    this->r = r; this->bg = bg; this->fg = fg; this->alpha = alpha;
    this->invert = invert;
    this->zstart = zstart; this->zstep = zstep;
  }
};

void dilate(volume8 &v, int r, float alpha, uint8 bg, uint8 fg, bool invert, uint8 uval, int num_threads) {
  v.pad(r, bg);

  vector<DilateThread *> threads(num_threads);
  for(int T = 0; T < num_threads; T++) 
    threads[T] = new DilateThread(v, r, bg, fg, alpha, invert, T, num_threads);
  for(int T = 0; T < num_threads; T++) threads[T]->start();
  for(int T = 0; T < num_threads; T++) threads[T]->wait();
  for(int T = 0; T < num_threads; T++) v.set(threads[T]->res, uval);
  for(int t = 0; t < num_threads; t++) delete threads[t];

  v.unpad(r);
}

struct ScaleTricubicThread : public QThread {
  volume8 *src, *dst; 
  int startScanLine, stepScanLine;  
  ScaleTricubicThread(int start, int step, volume8 *src, volume8 *dst) { 
    this->src = src; this->dst = dst;  startScanLine = start; stepScanLine = step;
  }
  void run() { dst->scale_tricubic(*src, startScanLine, stepScanLine); }
};

// Used to scale to z axis dimension.
volume8 *scaleToZ(volume8 *src, float sxy, int num_threads) {
  if(fabs(sxy - 1.0f) < 1e-5) return src; // no scaling necssary
  if(sxy < 1) {
    cout << "BLUR XY SIGMA = " << px2sigma(2.0f/sxy) << endl;
    blur(src, DimX, px2sigma(2.0f/sxy), num_threads);
    blur(src, DimY, px2sigma(2.0f/sxy), num_threads);
  }
  int width  = float(src->width) * sxy, height = float(src->height) * sxy;
  volume8 *dst = new volume8(width, height, src->depth);
  vector<ScaleTricubicThread *> threads;
  for(int t = 0; t < num_threads; t++) threads.push_back(new ScaleTricubicThread(t, num_threads, src, dst));
  for(int t = 0; t < num_threads; t++) threads[t]->start();
  for(int t = 0; t < num_threads; t++) threads[t]->wait();
  for(int t = 0; t < num_threads; t++) delete threads[t];
  delete src; return dst;
}

// Scale's z-axis to match x and y axes.
volume8 *scaleToXY(volume8 *src, float sz, int num_threads) {
  if(fabs(sz - 1.0f) < 1e-5) return src; // no scaling necssary
  cout << "BLUR Z SIGMA = " << px2sigma(2.0f/sz) << endl;
  if(sz < 1) blur(src, DimZ, px2sigma(2.0f/sz), num_threads);
  int depth = float(src->depth) * sz;
  volume8 *dst = new volume8(src->width, src->height, depth);
  vector<ScaleTricubicThread *> threads;
  for(int t = 0; t < num_threads; t++) threads.push_back(new ScaleTricubicThread(t, num_threads, src, dst));
  for(int t = 0; t < num_threads; t++) threads[t]->start();
  for(int t = 0; t < num_threads; t++) threads[t]->wait();
  for(int t = 0; t < num_threads; t++) delete threads[t];
  delete src;  
  return dst;
}

// Used to scale x and y axes so they're matched.
volume8 *scaleUpXY(volume8 *src, float sx, float sy, int num_threads) {
  if(fabs(sx - 1.0f) < 1e-5 && fabs(sy - 1.0f) < 1e-5) return src; // no scaling necssary
  volume8 *dst = new volume8(sx *src->width, sy * src->height, src->depth);
  cout << "scalling up (x,y)" << endl;
  cout << "(src) vx=" << src->width << " " << "sy=" << src->height << endl;
  cout << "(dst) vx=" << dst->width << " " << "sy=" << dst->height << endl;
  vector<ScaleTricubicThread *> threads;
  for(int t = 0; t < num_threads; t++) threads.push_back(new ScaleTricubicThread(t, num_threads, src, dst));
  for(int t = 0; t < num_threads; t++) threads[t]->start();
  for(int t = 0; t < num_threads; t++) threads[t]->wait();
  for(int t = 0; t < num_threads; t++) delete threads[t];
  delete src;  return dst;
}


volume8 *scaleXYZ(volume8 *src, float sxyz, int num_threads) {
  if(fabs(sxyz - 1.0f) < 1e-5) return src; // no scaling necssary
  if(sxyz < 1) {
    cout << "sxyz=" << sxyz << endl;
    blur(src, DimX, px2sigma(1.0f/sxyz), num_threads); blur(src, DimY, px2sigma(1.0f/sxyz), num_threads);
    blur(src, DimZ, px2sigma(1.0f/sxyz), num_threads);
  }
  int width  = float(src->width) * sxyz, height = float(src->height) * sxyz, depth = float(src->depth) * sxyz;
  volume8 *dst = new volume8(width, height, depth);

  cout << width << " " << height << " " << depth << endl;

  vector<ScaleTricubicThread *> threads;
  for(int t = 0; t < num_threads; t++) threads.push_back(new ScaleTricubicThread(t, num_threads, src, dst));
  for(int t = 0; t < num_threads; t++) threads[t]->start();
  for(int t = 0; t < num_threads; t++) threads[t]->wait();
  for(int t = 0; t < num_threads; t++) delete threads[t];
  delete src; return dst;
}



struct EDTthread1 : public QThread {
  volume32 &D; volume8 &I; uint8 fg; int zstart, zstep;
  EDTthread1(volume32 &D_, volume8 &I_, uint8 fg, int zstart, int zstep) : D(D_), I(I_) {
    this->fg = fg; this->zstart = zstart; this->zstep = zstep;
  }
  void run() { D.ComputeEDTstep1(I, fg, zstart, zstep); }
};
struct EDTthread2 : public QThread {
  volume32 &D;
  int xstart, xstep;
  EDTthread2(volume32 &D_, int xstart, int xstep) : D(D_), xstart(xstart), xstep(xstep) { }
  void run() { D.ComputeEDTstep2(xstart, xstep); }
};
// Multithreaded Euclidan distance transform (EDT). 
void compute_edt(volume32 &D, volume8 &I, uint8 fg, int num_threads) {
  vector<EDTthread1 *> threads1(num_threads);
  for(int t = 0; t < num_threads; t++) threads1[t] = new EDTthread1(D, I, fg, t, num_threads);
  for(int t = 0; t < num_threads; t++) threads1[t]->start();
  for(int t = 0; t < num_threads; t++) threads1[t]->wait();
  for(int t = 0; t < num_threads; t++) delete threads1[t];
  vector<EDTthread2 *> threads2(num_threads);
  for(int t = 0; t < num_threads; t++) threads2[t] = new EDTthread2(D, t, num_threads);
  for(int t = 0; t < num_threads; t++) threads2[t]->start();
  for(int t = 0; t < num_threads; t++) threads2[t]->wait();
  for(int t = 0; t < num_threads; t++) delete threads2[t];
}

// Applies rank filtering.
struct RankThread : public QThread {
  volume8 &dst, &paddedsrc;
  int r, startScanLine, stepScanLine;  
  float alpha;
  RankThread(volume8 &dst_, volume8 &paddedsrc_, int r_, float alpha_, int start, int step) : 
    dst(dst_), paddedsrc(paddedsrc_) { 
    r = r_;
    alpha = alpha_;
    startScanLine = start; stepScanLine = step; 
  }
  void run() { paddedsrc.rankfilt(dst, r, alpha, startScanLine, stepScanLine); }
};

void rankfilt(volume8 &v, int r, float alpha, int num_threads) {
  vector<RankThread *> threads(num_threads);
  volume8 padded(v.width + 2 * r, v.height + 2 * r, v.depth + 2 * r);
  padded.fill(0);
  // Copy binary volume into the padded volume.
  for(int z = 0; z < v.depth; z++) 
    for(int y = 0; y < v.height; y++) 
      for(int x = 0; x < v.width; x++) padded(x+r,y+r,z+r) = v(x,y,z);	  

  for(int T = 0; T < num_threads; T++) threads[T] = new RankThread(v, padded, r, alpha, T, num_threads);
  for(int t = 0; t < num_threads; t++) threads[t]->start();
  for(int t = 0; t < num_threads; t++) threads[t]->wait();
  for(int t = 0; t < num_threads; t++) delete threads[t];
}

// Collects points to thin.
struct ThinPtlThread : public QThread {
  vector<ivec3> ptl;
  volume8 &v;
  thindir_t thindir; thintype_t thintype;
  int zstart, zstep;
  ThinPtlThread(volume8 &v, thindir_t thindir, thintype_t thintype, int zstart, int zstep) : 
    v(v), thindir(thindir), thintype(thintype), zstart(zstart), zstep(zstep) { }
  void run() { v.thinptl(ptl, thindir, thintype, zstart, zstep); }
};

// Multithreaded thinning
void thin(volume8 &v, thintype_t thintype, int num_threads) {
  v.thin_init(); // Clears 1 voxel border for thinning.
  int modified = 0;
  do {
    modified = 0;
    // Iterate through each of 6 directions.
    for(int d = 0; d < 6; d++) {
      thindir_t thindir = (thindir_t)d;
      // First get point list for thinning.
      vector<ThinPtlThread *> threads1(num_threads);
      for(int t = 0; t < num_threads; t++) threads1[t] = new ThinPtlThread(v, thindir, thintype, t, num_threads);
      for(int t = 0; t < num_threads; t++) threads1[t]->start();
      for(int t = 0; t < num_threads; t++) threads1[t]->wait();

      // Single thread apply thinned points.
      for(int t = 0; t < num_threads; t++) modified += v.thin_apply(threads1[t]->ptl, thintype);
      for(int t = 0; t < num_threads; t++) delete threads1[t]; 
    }
  } while(modified > 0);
}


