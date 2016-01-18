#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector> 
#include <sstream>

#include <map>
#include <set>
#include <queue>

#include <tiffio.h>

#include "volume.hpp"
#include "util.hpp"
#include "marchingcubes.hpp"

namespace vol {
  // Constructs a sphereical structural element centered at (0,0,)
  // where the sphere centered at (cx,cy,cz) is removed.
  inline void sphere_strel(vector<ivec3> &strel, int r, int cx, int cy, int cz) {
    int rsq = r*r;
    int r2 = r + max(abs(cx), max(abs(cy), abs(cz)));
    for(int z = -r2; z <= r2; z++) 
      for(int y = -r2; y <= r2; y++) 
	for(int x = -r2; x <= r2; x++) {
	  // In sphere centered at (x,y,z).
	  if((x * x + y * y + z * z <= rsq) == true) {
	    // But, not in sphere centered at (cx,cy,cz).
	    int dx = x - cx, dy = y - cy, dz = z - cz;
	    if((dx * dx + dy * dy + dz * dz <= rsq) == false) {
	      strel.push_back(ivec3(x,y,z));
	    }
	  }
	}
  }

  // Copy a block of memory of length N from src to dst (vectorizable code).
  void fastcopy1(float *__restrict__  dst, uint8 *__restrict__ src, long N) { for(long i = 0; i < N; i++) dst[i] = src[i]; }
  void fastcopy2(uint8 *__restrict__  dst, float *__restrict__ src, long N) { for(long i = 0; i < N; i++) dst[i] = src[i]; }

  // Zero a bloc of memory (vectorizable code).
  void zero(float *x, long N) { for(long i = 0; i < N; i++) x[i] = 0; }
  // Fill block with given value.
  void fillval(float *x, long N, float s) { for(long i = 0; i < N; i++) x[i] = s; }

  // Perform convolution between f and x, storing product in acc. Note
  // that these vectors must be the same size |f| = |acc| = K.
  float convolve(float *__restrict__ f, float *__restrict__ x, float *__restrict__ acc, int K) {
    for(int k = 0; k < K; k++) acc[k] = f[k] * x[k];
    float sum = 0; for(int k = 0; k < K; k++) sum += acc[k];
    return sum;
  }

  // Uses libtiff to load a given frame from a hyperstack with the given slices per stack.
  bool volume16::load(string filename, int slices, int num_channels, int frame, int channel, int z_depth, int z_start) {
    TIFFSetWarningHandler(NULL); TIFFSetErrorHandler(NULL);
    TIFF *tif = TIFFOpen(filename.c_str(), "r");
    if(tif == NULL) return false;

    // Advance through image stack to given frame.
    for(int pos = 0; pos < (slices * num_channels) * frame; pos++) 
      if(!TIFFReadDirectory(tif)) { cerr << "pos=" << pos << " unable to advance to position" << endl; return false; }

    // Load image information. 
    uint16 bps = 0;
    if(z_depth > 0) depth = z_depth;
    else { depth = slices; z_start = 0; } // depth is the number of slices
    for(int z = z_start; z < depth; z++) {
      for(int c = 0; c < num_channels; c++) {
	if(c == channel) { // Read image from requested channel.
	  uint32 imageWidth;   TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imageWidth);
	  uint32 imageLength;  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imageLength);
	  uint16 bitsPerSample;TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitsPerSample); 
	  if((bitsPerSample == 8 || bitsPerSample == 16) == false) {
	    cerr << "invalid bits per sample" << endl; return false; 
	  }
	  is8bitOriginal = bitsPerSample == 8;
	  if(data == NULL) {
	    bps = bitsPerSample;
	    width = imageWidth; height = imageLength;
	    data = new uint16[width*height*depth];
	    setodim(width, height, depth); // set octree dimensions 
	  }
	  else {
	    if((int)imageWidth != width || (int)imageLength != height || bps != bitsPerSample) { 
	      cerr << "invalid image in stack" << endl; return false; 
	    }
	  }
	  // Read data branching on data type (floating point, gray8, gray16)
	  int offz  = z * width * height;
	  tdata_t buf = _TIFFmalloc(TIFFScanlineSize(tif));
	  for (uint32 row = 0; row < imageLength; row++) {
	    TIFFReadScanline(tif, buf, row);
	    int off = offz + row * width;
	    if(bitsPerSample == 8) { // Load a 8-bit gray scale file.
	      unsigned char *gray8 = (unsigned char*)buf;
	      for(int i = 0; i < width; i++) data[off + i] = gray8[i];
	    }
	    else { // Load a 16-bit gray scale file.
	      unsigned short *gray16 = (unsigned short *)buf;
	      for(int i = 0; i < width; i++) data[off + i] = gray16[i];
	    }
	  }
	  _TIFFfree(buf);
	}
	if(!TIFFReadDirectory(tif) && z != slices - 1) { cerr << "unable to advance to next slice" << endl; return false;}
      }
    }
    TIFFClose(tif); return true;
  }

  // Save function for 8-bit volumes.
  bool volume8::save(string filename) {
    TIFF *tif = TIFFOpen(filename.c_str(), "w");
    string raw_description;
    raw_description += "frame=1\n";
    raw_description += "slices=" + util::toString(depth) + "\n";
    raw_description += "channel=1\n";
    TIFFSetField(tif, TIFFTAG_IMAGEDESCRIPTION, raw_description.c_str()); 
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    for(int z = 0; z < depth; z++) {
      // Set only a minimal set of TIFF attributes, expand later.
      TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
      TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
      TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);

      tdata_t buf = _TIFFmalloc(width * sizeof(uint8));
      uint8 *gray8 = (uint8 *)buf;
      for(int y = 0; y < height; y++) {
	for(int x = 0; x < width; x++) gray8[x] = v(x,y,z);
	TIFFWriteScanline(tif, buf, y);
      }
      _TIFFfree(buf);
      TIFFWriteDirectory(tif);
    }
    TIFFClose(tif); return true;    
  }

  bool volume8::save_mha(string filename, float dx, float dy, float dz) {
    ostringstream header;

    header << "ObjectType = Image" << endl;
    header << "NDims = 3" << endl;
    header << "BinaryData = True" << endl;
    header << "BinaryDataByteOrderMSB = False" << endl;
    header << "CompressedData = False" << endl;
    header << "TransformMatrix = 1 0 0 0 1 0 0 0 1"  << endl;
    header << "Offset = 0 0 0" << endl;
    header << "CenterOfRotation = 0 0 0" << endl;
    header << "AnatomicalOrientation = RAI" << endl;
    header << "ElementSpacing = " << dx << " " << dy << " " << dz << endl;
    header << "DimSize = " << width << " " << height << " " << depth << endl;
    header << "ElementType = MET_UCHAR" << endl;
    header << "ElementDataFile = LOCAL" << endl;

    string header_str = header.str();
    size_t header_len = header_str.length();

    ofstream mha(filename.c_str(), ofstream::binary);
    mha.write((char*)&header_str[0], header_len);
    mha.write((char*)data, width*height*depth);
    mha.close();
    return true;
  }

  // Compute a linear Gaussian filter with the given sigma. fdim is
  // set to the filter dimension which is 2 * dim + 1.
  float *Gaussian(float sigma, int &dim, int &fdim) {
    dim = (int)ceil(2 * sigma); 
    if(dim <= 1) dim = 1;
    fdim = 2*dim + 1; 
    float *filter = new float[fdim];

    // Compute the filter.
    float sigmasq = sigma * sigma, Z = 0;
    for(int i = -dim; i <= dim; i++) {
      float Fi = exp( -0.5 * i * i / sigmasq);
      Z += Fi; 
      filter[i + dim] = Fi;
    }
    // Normalize the Gaussian filter.
    for(int i = 0; i < fdim; i++) filter[i] /= Z;
    return filter;
  }

  void volume8::GaussX(float sigmaX, int ystart, int ystep) {
    int dim, fdim; float *filter = Gaussian(sigmaX, dim, fdim); 

    float *buf = new float[fdim];
    float *src = new float[width + 2 * dim]; 
    float *dst = new float[width];
    for(int z = 0; z < depth; z++) {
      int offz = z * width * height;
      for(int y = ystart; y < height; y += ystep) { 
	int offy = offz + y * width;
	// Copy edge pixels.
	// v(x = 0, y, z)
	fillval(src, dim, *(data + offy));
	// v(x=width-1, y, z)
	fillval(src + width + dim, dim, *(data + offy + width - 1));
	// Copy actual data.
	fastcopy1(src + dim, data + offy, width);
	for(int x = 0; x < width; x++)  dst[x] = convolve(filter, src + x, buf, fdim);
	fastcopy2(data + offy, dst, width);
      }
    }

    delete[] src; delete[] dst; delete[] buf; delete[] filter;
  }

  void volume8::GaussY(float sigmaY, int xstart, int xstep) {
    int dim, fdim; float *filter = Gaussian(sigmaY, dim, fdim);

    float *buf = new float[fdim];
    // src is bigger to allow for over run on the sides.
    float *src = new float[height + 2 * dim]; fillval(src, height + 2 * dim, 0.0);
    float *dst = new float[height];

    for(int z = 0; z < depth; z++) {
      int offz = z * width * height;
      for(int x = xstart; x < width; x += xstep) {
	// Copy y-axis scan line.
	for(int y = 0; y < height; y++) src[y + dim] = data[offz + y * width + x];

	// Copy edge pixels to border.
	// v(x,y=0,z)
	float d0 = data[offz + x]; for(int s = 0; s < dim; s++) src[s] = d0;
	// v(x,y=height-1,z)
	float d1 = data[offz + (height-1) * width + x]; for(int s = 0;  s < dim; s++) src[height+dim+s] = d1;
	// Perform convolution.
	for(int y = 0; y < height; y++) dst[y] = convolve(filter, src + y, buf, fdim);
	// Copy buffer back.
	for(int y = 0; y < height; y++) data[offz + y * width + x] = dst[y];
      }
    }    
    delete[] src; delete[] dst; delete[] buf; delete[] filter;
  }

  void volume8::GaussZ(float sigmaZ, int xstart, int xstep) {
    int dim, fdim; float *filter = Gaussian(sigmaZ, dim, fdim);

    float *buf = new float[fdim];
    // src is bigger to allow for over run on the sides.
    float *src = new float[depth + 2 * dim]; fillval(src, depth + 2 * dim, 0.0);
    float *dst = new float[depth];

    int wh = width * height;
    for(int y = 0; y < height; y++) {
      int offy = y * width;
      for(int x = xstart; x < width; x += xstep) {
	// v(x,y,z=0)
	float d0 = data[offy+x]; for(int s = 0; s < dim; s++) src[s] = d0;
	// v(x,y,z=depth-1)
	float d1 = data[wh*(depth-1)+offy+x]; for(int s = 0; s < dim; s++) src[depth+dim+s] = d1;
	// Copy z-axis scan line.
	for(int z = 0; z < depth; z++) src[z + dim] = data[ wh * z + offy + x];
	// Perform convolution.
	for(int z = 0; z < depth; z++) dst[z] = convolve(filter, src + z, buf, fdim);
	// Overwrite z-axis scan

	for(int z = 0; z < depth; z++) data[wh * z + offy + x] = dst[z];
      }
    }
    delete[] src; delete[] dst; delete[] buf; delete[] filter;
  }

  void volume8::dfs(volumeT<char> &labels, ivec3 start, vector<ivec3> &comp, unsigned char fg, bool use_nhood26) {
    vector<ivec3> L;  L.reserve(400); L.push_back(start); 
    labels(start) = 1; // 1 = in queue
    vector<ivec3> adj; adj.reserve(26);
    while(L.empty() == false) {
      // Get current voxel from queue.
      ivec3 cur = L.back(); L.pop_back();
      labels(cur) = 2; // mark as finished
      comp.push_back(cur);
      adj.clear(); 
      if(use_nhood26) nhood26(adj, cur); else nhood6(adj, cur);       // Get adjacent voxels.
      for(size_t j = 0; j < adj.size(); j++) {
	// Has it been visited and is it the fg value? 
	if(labels(adj[j]) == 0 && v(adj[j]) == fg) {
	  labels(adj[j]) = 1; // label 1 to indicate in stack
	  L.push_back(adj[j]);  // add to stack
	}
      }
    }
  }

  // Returns connected components.
  void volume8::components(vector< vector<ivec3> > &comps, int minvxls, unsigned char fg, bool use_nhood26) {
    volumeT<char> labels(width, height, depth);
    labels.fill(0);
    vector<ivec3> comp;
    for(int z = 0; z < depth; z++) 
      for(int y = 0; y < height; y++) 
	for(int x = 0; x < width; x++) {
	  if(labels(x,y,z) != 0) continue; // already in a component, skip
	  if(v(x,y,z) != fg) continue; // not fg value
	  comp.clear(); 
	  dfs(labels, ivec3(x, y, z), comp, fg, use_nhood26);
	  if(minvxls <= (int)comp.size()) comps.push_back(comp); 
	}
  }

  // Gets foreground voxels adjacent to background.
  void volume8::outline_vox(vector<ivec3> &ovox, uint8 fg, bool use_nhood26) {
    vector<ivec3> hood(26);
    for(int z = 0; z < depth; z++) 
      for(int y = 0; y < height; y++) 
	for(int x = 0; x < width; x++) {
	  if(v(x,y,z) != fg) continue;
	  // Get neighborhood.
	  if(use_nhood26) nhood26(hood, ivec3(x,y,z)); else nhood6(hood, ivec3(x,y,z));
	  for(size_t h = 0; h < hood.size(); h++) {
	    if(v(hood[h]) != fg) { ovox.push_back(ivec3(x,y,z)); break; }
	  }
	}
  }

  void volume8::triangulate(vector<vec3> &vtable, vector<face> &ftable, int mcdim, uint8 fg, uint8 bg) {
    float fgf = fg, bgf = bg;
    // This is absolutely critical for the shift/alpha code to work properly below.
    // Want the isosurface exactly half-way in between FG and BG values
    float isolevel = (fgf - bgf) / 2.0; 
    pad(mcdim, 0);
    // TODO: Move MarchingCubes to volumeT<T>!!!
    marching::MarchingCubes(isolevel, *this, vtable, ftable, mcdim);    
    unpad(mcdim);

    // This shift assures that voxelize in geom::mesh returns almost
    // the same same binary volume.
    float shift = float(mcdim) / 2.0f; // The shift gets rid of "left-most" cube. 
    vec3 alpha;
    alpha[0] = (float(width )) / float(width  + 1); // alpha scales to compensate for the second cube added at the end
    alpha[1] = (float(height)) / float(height + 1); // during pad() above
    alpha[2] = (float(depth )) / float(depth  + 1);
    for(size_t i = 0; i < vtable.size(); i++) { vtable[i] -= shift; vtable[i] *= alpha; }
  }

  inline void volume8::addstrel(vector<int> &hist, vector<ivec3> &s, int x, int y, int z) {
    for(size_t i = 0; i < s.size(); i++) 
      hist[v(x + s[i].x(), y + s[i].y(), z + s[i].z())]++;
  }

  inline void volume8::substrel(vector<int> &hist, vector<ivec3> &s, int x, int y, int z) {
    for(size_t i = 0; i < s.size(); i++) hist[v(x + s[i].x(), y + s[i].y(), z + s[i].z())]--;
  }

  void volume8::dilate(vector<ivec3> &update, 
		       int r, uint8 fg, uint8 bg, float alpha, bool invert,
		       int zstart, int zstep) {
    vector<ivec3> s, sL, sR, sU, sD, sT, sB; 
    sphere_strel(s, r); // entire sphere centered at (0,0,0)
    int thresh = s.size() * alpha;
    // Prepare structural elements for updating fg value count.
    sphere_strel(sL, r,  1,  0, 0);     // remove sphere to right
    sphere_strel(sR, r, -1,  0, 0);     // remove sphere to left
    sphere_strel(sU, r,  0,  1, 0);     // remove sphere from below
    sphere_strel(sD, r,  0,  -1, 0);    // remove sphere from above
    sphere_strel(sT, r,  0, 0, zstep);  // remove sphere from top
    sphere_strel(sB, r,  0, 0, -zstep); // remove sphere from bottom

    int x = r, y = r, z = r + zstart, xstep = 1, ystep = 1; 
    vector<int> hist(256, 0);
    addstrel(hist, s, x, y, z);
    int d = depth - r, w = width - r, h = height - r;
    while(z < d) {
      while(r <= y && y < h) {
	while(r <= x && x < w) {
	  if(v(x,y,z) == bg) {
	    if(invert == false) { if(hist[fg] > thresh)  update.push_back(ivec3(x,y,z)); }
	    else { if(hist[fg] <= thresh) update.push_back(ivec3(x,y,z)); }
	  }
	  x += xstep;
	  if((xstep > 0 && x < w) || (xstep < 0 && x >= r)) {
	    // Depending on direction of scan, pick correct STREL to
	    // use for updating fg counts.
	    vector<ivec3> &s1 = xstep > 0 ? sL : sR;
	    substrel(hist, s1, x-xstep,y,z);
	    vector<ivec3> &s2 = xstep > 0 ? sR : sL;
	    addstrel(hist, s2, x,y,z);
	  }
	}
	xstep = -xstep; x += xstep; // flip direction of x-scan
	y += ystep;
	if((ystep > 0 && y < h) || (ystep < 0 && y >= r)) {
	  // Depending on the direction of the scan pick the corresponding
	  // strel for updating the y position.
	  vector<ivec3> &s1 = ystep > 0 ? sU : sD;
	  substrel(hist, s1, x,y-ystep,z);
	  vector<ivec3> &s2 = ystep > 0 ? sD : sU;
	  addstrel(hist, s2, x,y,z);
	}
      }
      ystep = -ystep; y += ystep; // flip direction of y-scan.
      z += zstep;
      if( z < d ) {
	// The z-scan always goes deeper into the image volume.
	substrel(hist, sT, x,y,z-zstep);
	addstrel(hist, sB, x,y,z);
      }
    }      
  }

  bool volume32::load_mha_float(const string &filename) {
    ifstream mha(filename.c_str(), ios::binary | ios::ate);
    mha.seekg(-(sizeof(float) * width * height * depth), ios::end);
    float *buf = new float[width * height * depth];
    mha.read((char*)buf, sizeof(float) * width * height * depth);
    for(int i = 0; i < width*height*depth; i++) {
      data[i] = buf[i];
    }
    delete[] buf;
    return true;
  }

  inline bool volume32::RemoveEDT(int d2u, int d2v, int d2w, int ud, int vd, int wd) {
    int a = vd - ud, b = wd - vd, c = wd - ud;
    return c * d2v - b * d2u - a * d2w - a * b * c > 0;
  }
  void volume32::VoronoiEDT(int d, ivec3 j) {
    int maxnd = max(width, max(height, depth));
    vector<int> h(maxnd+1), g(maxnd+1);
    int l = 0, nd = 0;
    if(d == 1) nd = width; else if(d == 2) nd = height; else if(d == 3) nd = depth;
    ivec3 xi = j;
    for(int i = 0; i < nd; i++) {
      xi[d-1] = i;
      int fi = D(xi);
      if( fi != numeric_limits<int>::max() ) {
	if(l < 2) {  l = l+1; g[l] = fi; h[l] = i; }
	else {
	  while( l >= 2 && RemoveEDT(g[l-1], g[l], fi, h[l-1], h[l], i) ) l = l - 1;
	  l = l + 1; g[l] = fi; h[l] = i;
	}
      }
    }
    int ns = l; if(ns == 0) return;
    l = 1;
    for(int i = 0; i < nd; i++) {
      xi[d-1] = i;
      while( l < ns && g[l] + sqr(h[l] - i) > g[l+1] + sqr(h[l+1] - i) ) l = l + 1;
      D(xi)  = g[l] + sqr(h[l] - i);
    }
  }
  void volume32::ComputeEDT(volume8 &I, uint8 fg, int d, ivec3 j) {
    if(d == 1) {
      for(int x = 0; x < width; x++) {
	if(I(x, j[1], j[2]) == fg) D(x, j[1], j[2]) = 0;
	else D(x, j[1], j[2]) = numeric_limits<int>::max();
      }
    }
    else {
      if(d == 2) { for(int y = 0; y < height; y++) ComputeEDT(I, fg, d - 1, ivec3(-1, y, j[2])); }
      else if(d == 3) { for(int z = 0; z < depth;  z++) ComputeEDT(I, fg, d - 1, ivec3(-1,-1,z) ); }
    }
    if(d == 1) VoronoiEDT(d, ivec3(-1, j[1], j[2]));
    else if(d == 2) {
      for(int x = 0; x < width; x++) VoronoiEDT(d, ivec3(x, -1, j[2]));
    }
    else if(d == 3) {
      for(int x = 0; x < width; x++)
	for(int y = 0; y < height; y++) VoronoiEDT(d, ivec3(x,y,-1));
    }
  }

  // Fill volume with 0 val, and seed regions in volume.
  void volume32::seed(vector< vector<ivec3> > &seeds) {
    fill(0);
    for(int i = 0; i < (int)seeds.size(); i++) {
      int label = i + 1;
      vector<ivec3> &seed = seeds[i];
      for(size_t s = 0; s < seed.size(); s++) v(seed[s][0], seed[s][1], seed[s][2]) = label;
    }
  }

  // Adapted from "The Watershed Transform in ITK - discussions and
  // new developments" by Richard Beare and Gaeten Lehmann.  This is
  // an algorithm attributed to Beucher.
  // Convention for labeled image.
  // O(x,y,z) = 0, unlabeled; O(x,y,z) = 1... S, seed regions
  // O(x,y,z) <= -1, do not touch this region.
  void volume32::watershed(volume32 &I, bool use_nhood26) {
    map<int, queue<ivec3> > PQ;
    vector<ivec3> N; // neighborhood positions
    
    for(int z = 0; z < depth; z++)
      for(int y = 0; y < height; y++)
	for(int x = 0; x < width; x++) {
	  ivec3 p(x,y,z);
	  if(O(p) > 0) {
	    // Put all marker/seed pixels with background neighbors into
	    // the priority queue.
	    bool bgNeighbor = false;
	    N.clear();
	    if(use_nhood26) nhood26(N,p); else nhood6(N,p);
	    for(size_t i = 0; i < N.size(); i++) {
	      ivec3 &q = N[i];
	      if(O(q) == 0) bgNeighbor = true;
	    }
	    if(bgNeighbor) PQ[I(p)].push(p); // Have separate queues for each intensity level.
	  }
	}

    while(!PQ.empty()) {
      // Get the queue for the lowest intensity level.
      int vimg = PQ.begin()->first; 
      queue<ivec3> curQ = PQ.begin()->second;
      PQ.erase( PQ.begin() );
      while(curQ.empty() == false) { // Empty this queue.
	ivec3 p = curQ.front(); curQ.pop();
	// Propagate the label of p to all unlabeled neighbors of p.
	N.clear(); 
	if(use_nhood26) nhood26(N, p); else nhood6(N, p); 
	for(size_t i = 0; i < N.size(); i++) {
	  ivec3 &q = N[i];  // Check to see if q is an unlabeled neighbor.
	  if(O(q) == 0) {  // Propagate the label of p to unlabeled neighbors.
	    O(q) = O(p); 
	    // Is the intensity at the propagated position less than
	    // (or equal to) the current queue intensity? Add to
	    // current queue. This deals "reasonably" with plateaus in
	    // the image.
	    if ( I(q) <= vimg ) curQ.push(q); 
	    else PQ[I(q)].push(q); // Otherwise add to corresponding intensity queue.
	  }
	}
      } 
    }

  }

  void volume32::components(vector< vector<ivec3> > &final, bool clear_empty) {
    long whd = width * height * depth;
    int num_comp = -1;
    vector< vector<ivec3> > comps;
    for(long i = 0; i < whd; i++) num_comp = max(num_comp, data[i]);
    comps.resize(num_comp);
    for(int z = 0; z < depth; z++) 
      for(int y = 0; y < height; y++)
	for(int x = 0; x < width; x++) {
	  if(v(x,y,z) > 0) comps[v(x,y,z) - 1].push_back(ivec3(x,y,z));
	}
    if(clear_empty) {
      for(size_t c = 0; c < comps.size(); c++) { if(comps[c].size() > 0) final.push_back(comps[c]); }
    }
    else final = comps;
  }
  
  // Fractional Hausdorff distance. Assumes volume has squared EDT values.
  float volume32::frac_hausdorff(const vector<ivec3> &vox, float alpha) {
    if(vox.size() == 0) { return -1; }
    vector<int> EDT_vals;
    EDT_vals.reserve(vox.size());
    // Get values from EDT.
    for(size_t i = 0; i < vox.size(); i++) {
      int x = vox[i][0], y = vox[i][1], z = vox[i][2];
      // Make sure voxel point is in range.
      if( 0 <= x && x < width && 0 <= y && y < height && 0 <= z && z < depth ) 
	EDT_vals.push_back(v(x,y,z));
    }   
    // If the entire template is out of range, return max distance.
    if(EDT_vals.size() == 0) return -1;

    // Sort and return fractional distance value.
    sort(EDT_vals.begin(), EDT_vals.end());
    float sum = 0, N = EDT_vals.size() * alpha; // also trying trimmed mean
    for(size_t i = 0; i < EDT_vals.size() * alpha; i++) sum += EDT_vals[i] / N;
    return sum;
  }

  // For tracking, you need some way of resolving targets based on distance.
  // Tracker will "jump" onto neighbor if best match.
  float volume32::frac_hausdorff_match(ivec3 &min_shift, vector<ivec3> &vox, int d, float alpha) {
    int dsq = d * d;
    float min_dist = numeric_limits<float>::max();
    bool found_min = false;
    vector<ivec3> vox_shift = vox;
    for(int dx = -d; dx <= d; dx++) 
      for(int dy = -d; dy <= d; dy++)
	for(int dz = -d; dz <= d; dz++) 
	  // Limit shifts to a spherical region.
	  if(dx * dx + dy * dy + dz * dz <= dsq) {
	    // Allow shift and get fractional hausdorf distance.
	    ivec3 shift(dx,dy,dz);
	    for(size_t i = 0; i < vox.size(); i++) vox_shift[i] = vox[i] + shift;
	    float dist = frac_hausdorff(vox_shift, alpha);
	    if(dist >= 0 && dist < min_dist) {
	      min_dist = dist;
	      min_shift = shift;
	      found_min = true;
	    }
	  }
    if(found_min == false) return -1;  else return min_dist;
  }

  void volume32::max_overlap(int &max_label, int &max_cnt, vector<ivec3> &vox) {
    // Get the labels and their overlap.
    vector<int> labels, label_cnts;
    for(size_t i = 0; i < vox.size(); i++) {
      if(vox[i][0] < 0 || vox[i][1] < 0 || vox[i][2] < 0) continue;
      if(vox[i][0] >= width || vox[i][1] >= height || vox[i][2] >= depth) continue;

      int cur_label = v(vox[i]);
      //  if(cur_label == 0) continue; // Make this a parameter?
      bool found_label = false;
      for(size_t l = 0; l < labels.size(); l++) { 
	// Find label, update count.
	if(labels[l] == cur_label) { label_cnts[l]++; found_label = true; break; } 
      }
      if(!found_label) { label_cnts.push_back(1);  labels.push_back(cur_label); }
    }
    max_cnt = 0; max_label = -1;
    for(size_t l = 0; l < label_cnts.size(); l++) { 
      if(label_cnts[l] > max_cnt) { max_cnt = label_cnts[l]; max_label = labels[l]; }
    }
  }

  void volume32::frac_hausdorff_match(vector<Hmatch_t> &matches, vector<ivec3> &vox, vector<ivec3> &voxfill, int d, float alpha, volume32 &labels) {
    int dsq = d * d;
    vector<ivec3> vox_shift = vox, voxfill_shift = voxfill;
    for(int dx = -d; dx <= d; dx++) 
      for(int dy = -d; dy <= d; dy++)
	for(int dz = -d; dz <= d; dz++) 
	  // Limit shifts to a spherical region.
	  if(dx * dx + dy * dy + dz * dz <= dsq) {
	    // Apply shift.
	    ivec3 shift(dx,dy,dz);
	    for(size_t i = 0; i < vox.size(); i++) vox_shift[i] = vox[i] + shift;
	    for(size_t i = 0; i < voxfill_shift.size(); i++) voxfill_shift[i] = voxfill[i] + shift;

	    float dist = frac_hausdorff(vox_shift, alpha);
	    int cur_max_label, cur_max_cnt;
	    labels.max_overlap(cur_max_label, cur_max_cnt, voxfill_shift);

	    bool found_label = false;
	    for(size_t l = 0; l < matches.size(); l++) {
	      if(matches[l].max_label == cur_max_label) {
		found_label = true;
		if(dist < matches[l].dist) { 
		  matches[l].shift = shift; 
		  matches[l].dist = dist; 
		  matches[l].max_overlap = cur_max_cnt;
		}
	      }
	    }
	    if(!found_label) { // New max label found.
	      Hmatch_t newmatch;
	      newmatch.shift = shift; 
	      newmatch.dist = dist;
	      newmatch.max_label = cur_max_label;
	      newmatch.max_overlap = cur_max_cnt;
	      matches.push_back(newmatch);
	    }
	  }
    // sort matches on distance, smallest to largest
    sort(matches.begin(), matches.end());
  }

  void volume32::hystersis_threshold(volume8 &hout, int thresh_low, int thresh_high) {
    if(thresh_low > thresh_high) { cerr << "hystersis_threshold(): thresh_low > thresh_high" << endl; exit(1); }
    volume8 high_thresh(width, height, depth);
    hout.fill(0);
    high_thresh.fill(0);
    for(int z = 0; z < depth; z++) 
      for(int y = 0; y < height; y++) 
	for(int x = 0; x < width; x++) {
	  if(v(x,y,z) >= thresh_low) hout(x,y,z) = 255;
	  if(v(x,y,z) >= thresh_high) high_thresh(x,y,z) = 255;
	}
    binary_reconstruct(hout, high_thresh);
  }

  // Code for loading in MATLAB:
  // alpha = XXX
  // w = XXX
  // h = XXX
  // d = XXX
  // fid = fopen("labelfile.binary")
  // label_matrix = fread(fid, w * h * d, "*int32")
  // labels3d = reshape(label_matrix, [w h d])
  // fclose(fid)
  void volume32::save_raw(const string &filename) {
    ofstream binary(filename.c_str(), ios::out | ios::binary);
    long datasz = width * height * depth;
    binary.write((char*)data, datasz * sizeof(int) );
  }

  // TODO: Precompute each of these cube 2^26 configurations and the output, store in a table.
  struct cube3x3_18 {
    uint8 v[3][3][3];
    ivec3 idx2pos[18];
    int pos2idx[3][3][3];
    int N6[18][6];
    vector<int> L; // stack for DFS
    vector<int> N6labels;
    uint8 & operator () (int x, int y, int z) { return v[x+1][y+1][z+1]; } 
    cube3x3_18() {
      int idx = 0;
      for(int z = 0; z < 3; z++)
	for(int x = 0; x < 3; x++)
	  for(int y = 0; y < 3; y++) {
	    // Don't use center voxel.
	    if(x == 1 && y == 1 && z == 1) pos2idx[x][y][z] = -1;
	    else {
	      // Only keep the N18 set.
	      if(z == 0 || z == 2) {
		if(x == 0 && y == 0) { pos2idx[x][y][z] = -1; continue; }
		if(x == 2 && y == 2) { pos2idx[x][y][z] = -1; continue; }
		if(x == 0 && y == 2) { pos2idx[x][y][z] = -1; continue; }
		if(x == 2 && y == 0) { pos2idx[x][y][z] = -1; continue; }
	      }
	      idx2pos[idx] = ivec3(x,y,z);
	      pos2idx[x][y][z] = idx;
	      idx++;
	    }
	  }

      // Clear neighbor table.
      for(int i = 0; i < 18; i++) for(int j = 0; j < 6; j++) N6[i][j] = -1;

      // Get 6-connected neighbors in sub cube.
      for(int i = 0; i < 18; i++) {
	int x = idx2pos[i][0], y = idx2pos[i][1], z = idx2pos[i][2];
	for(int dx = -1; dx <= 1; dx++)
	  for(int dy = -1; dy <= 1; dy++)
	    for(int dz = -1; dz <= 1; dz++) {
	      if(dx == 0 && dy == 0 && dz == 0) continue; // Skip center.
	      // Get only 6-connected neighbors
	      if((dx == 0 && dy == 0) || (dy == 0 && dz == 0) || (dx == 0 && dz == 0)) {
		int nx = x + dx, ny = y + dy, nz = z + dz;
		if( 0 <= nx && nx <= 2 &&  // Avoid outside of the cube.
		    0 <= ny && ny <= 2 && 
		    0 <= nz && nz <= 2 ) {
		  // Stay in 18-cube region.
		  if(nz == 0 || nz == 2) {
		    if(nx == 0 && ny == 0) continue; if(nx == 2 && ny == 2) continue;
		    if(nx == 0 && ny == 2) continue; if(nx == 2 && ny == 0) continue;
		  }
		  if( !(nx == 1 && ny == 1 && nz == 1) ) { // Skip center neighbor.
		    int neighidx = pos2idx[nx][ny][nz];
		    for(int j = 0; j < 6; j++) { if(N6[i][j] == -1) { N6[i][j] = neighidx; break; } }
		  }
		}
	      }
	    }	
      }
    }

    void dfs6(int labels[18], int lval, int start, uint8 bg = 0) {
      L.clear();
      L.push_back(start);
      while(L.empty() == false) {
	int cur = L.back(); L.pop_back();
	labels[cur] = lval;
	int j = 0;
	while( j < 6 && N6[cur][j] != -1 ) {
	  int nidx = N6[cur][j];
	  int x = idx2pos[nidx][0], y = idx2pos[nidx][1], z = idx2pos[nidx][2];
	  if(labels[nidx] == 0 && v[x][y][z] == bg) {
	    labels[nidx] = 1;
	    L.push_back(nidx);
	  }
	  j++;
	}
      }
    }
    int components18(uint8 bg = 0) {
      int labels[18]; for(int i = 0; i < 18; i++) labels[i] = 0; 
      int num_comps = 0;
      for(int idx = 0; idx < 18; idx++) {
	int x = idx2pos[idx][0], y = idx2pos[idx][1], z = idx2pos[idx][2];
	if(labels[idx] != 0) continue; // Already labeld vertex.
	if(v[x][y][z] != bg) continue;
	dfs6(labels, num_comps + 2, idx, bg);
	num_comps++;
      }
      N6labels.clear();
      for(int z = 0; z < 3; z++)
	for(int x = 0; x < 3; x++)
	  for(int y = 0; y < 3; y++) {
	    if(v[x][y][z] == bg) {
	      if((x == 1 && y == 1) || 
		 (y == 1 && z == 1) || 
		 (x == 1 && z == 1)) {
		int lcur = labels[pos2idx[x][y][z]];
		if(lcur == 0) continue; 
		bool found = false;
		for(size_t j = 0; j < N6labels.size(); j++) {
		  if(N6labels[j] == lcur) { found = true; break; }
		}
		if(!found) N6labels.push_back(lcur);
	      }
	    }
	  }
      return N6labels.size();
    }
    bool is_cond_4_satisfied(uint8 bg = 0) { return components18(bg) == 1; }
  };

  // TODO: Precompute 2^26 cube values and store in a table.
  struct cube3x3 {
    uint8 v[3][3][3];
    ivec3 idx2pos[26]; // 
    int pos2idx[3][3][3]; // cube position to index (1,1,1) == -1
    // N26 is a mapping for a cube index to a set of neighbor cube
    // indexes. Each are lists that end with -1.
    int N26[26][26];
    vector<int> L; // stack for DFS
    uint8 & operator () (int x, int y, int z) { return v[x+1][y+1][z+1]; } 
    cube3x3() {
      // Create forward and reverse mapping from position to an index in 3x3x3 cube.
      int idx = 0;
      for(int z = 0; z < 3; z++) 
	for(int x = 0; x < 3; x++) 
	  for(int y = 0; y < 3; y++) {
	    if(x == 1 && y == 1 && z == 1) pos2idx[x][y][z] = -1;
	    else {
	      idx2pos[idx] = ivec3(x,y,z); // index to position in 3-cube
	      pos2idx[x][y][z] = idx; // 3-cube position to index
	      idx++;
	    }
	  }

      for(int i = 0; i < 26; i++) for(int j = 0; j < 26; j++) N26[i][j] = -1;

      // Get valid 26-neighbors with center voxel removed.
      for(int i = 0; i <  26; i++) {
	int x = idx2pos[i][0], y = idx2pos[i][1], z = idx2pos[i][2];
	// Get 26-neighbors of an (x,y,z) position in a 3x3x3 cube.
	for(int dx = -1; dx <= 1; dx++)
	  for(int dy = -1; dy <= 1; dy++)
	    for(int dz = -1; dz <=1; dz++) {
	      if(dx == 0 && dy == 0 && dz == 0) continue; // Skip no-neighbor case.
	      // Only store neighbors within the 3x3x3 cube.
	      int nx = x + dx, ny = y + dy, nz = z + dz;
	      if( 0 <= nx && nx <= 2 && 
		  0 <= ny && ny <= 2 && 
		  0 <= nz && nz <= 2 ) {
		if( !(nx == 1 && ny == 1 && nz == 1) ) { // Skip center voxel.
		  int neighidx = pos2idx[nx][ny][nz];
		  for(int j = 0; j < 26; j++) { if( N26[i][j] == -1 ) { N26[i][j] = neighidx; break; } }
		}
	      }
	    }
      }
    }
    void dfs26(int labels[26], int lval, int start, uint8 fg = 255) {
      L.clear();
      L.push_back(start);
      while(L.empty() == false) {
	int cur = L.back(); L.pop_back();
	labels[cur] = lval;
	int j = 0;
	while( j < 26 && N26[cur][j] != -1 ) {
	  int nidx = N26[cur][j];
	  int x = idx2pos[nidx][0], y = idx2pos[nidx][1], z = idx2pos[nidx][2];
	  if(labels[nidx] == 0 && v[x][y][z] == fg) {
	    labels[nidx] = 1;
	    L.push_back(N26[cur][j]);
	  }
	  j++;
	}
      }
    }
    
    int components26(uint8 fg = 255) {
      int labels[26]; for(int i = 0; i < 26; i++) labels[i] = 0; 
      int num_comps = 0;
      for(int idx = 0; idx < 26; idx++) {
	if(labels[idx] != 0) continue; // Already labeled vertex.
	int x = idx2pos[idx][0],
	    y = idx2pos[idx][1],
	    z = idx2pos[idx][2];
	if(x == 1 && y == 1 && z == 1) continue; // skip center voxel
	if(v[x][y][z] != fg) continue;
	dfs26(labels, num_comps + 2, idx, fg);
	num_comps++;
      }
      return num_comps;
    }

    bool is_cond_2_satisfied(uint8 fg = 255) { return components26(fg) == 1; }
  };

  // Y(x,y,z) == 255, checks to see if point is == 0 in given direction specified
  // in the dir_t enum.
  bool volume8::is_border_point(thindir_t direction, int x, int y, int z) {
    switch(direction) {
    case thinU: return v(x, y, z - 1) == 0; case thinD: return v(x, y, z + 1) == 0;
    case thinN: return v(x, y - 1, z) == 0; case thinS: return v(x, y + 1, z) == 0;
    case thinE: return v(x - 1, y, z) == 0; case thinW: return v(x + 1, y, z) == 0;
    }
    return false;
  }

  bool volume8::is_endpoint(int x, int y, int z, uint8 fg) {
    int sum = 0;
    for(int dz = -1; dz <= 1; dz++) 
      for(int dy = -1; dy <= 1; dy++) 
	for(int dx = -1; dx <= 1; dx++) {
	  if(!(dx == 0 && dy == 0 && dz == 0)) { // Skip center.
	    if(v(x + dx, y + dy, z + dz) == fg) sum++;
	  }
	}
    if(sum > 1) return false; else return true;
  }
  // Used to 
  bool volume8::is_surface_endpoint(int x, int y, int z, uint8 fg) {
    // N6 contains one set of opposite white (background points)
    return 
      ( v(x-1, y, z) != fg && v(x+1, y, z) != fg ) ||
      ( v(x, y+1, z) != fg && v(x, y-1, z) != fg ) ||
      ( v(x, y, z-1) != fg && v(x, y, z+1) != fg );
  }
  bool volume8::is_finished(thintype_t thintype, int x, int y, int z, uint8 fg) {
    switch(thintype) {
    case thinSURFACE: return is_surface_endpoint(x, y, z, fg); 
    case thinLINE: return is_endpoint(x, y, z, fg); 
    default: cout << "is_finished()" << endl; exit(1);
    }
  }

  void volume8::thinptl(vector<ivec3> &ptl, thindir_t direction, thintype_t thintype, int zstart, int zstep) {
    uint8 fg = 255;
    cube3x3 cube26; cube3x3_18 cube6;
    for(int z = zstart+1; z < depth - 1; z += zstep) {
      for(int y = 1; y < height - 1; y++) 
	for(int x = 1; x < width - 1; x++)
	  if(v(x,y,z) == fg) { // Only need to process forground voxels.
	    if(is_border_point(direction, x, y, z)) {
	      if(!is_finished(thintype, x,y,z,fg)) {
		// Fill cubes for condition check if thinning not finished.
		for(int dz = -1; dz <= 1; dz++) 
		  for(int dy = -1; dy <= 1; dy++)
		    for(int dx = -1; dx <= 1; dx++) {
		      cube26(dx,dy,dz) = v(x+dx,y+dy,z+dz); cube6(dx,dy,dz) =  v(x+dx,y+dy,z+dz);
		    }
		// Check thinning conditions.
		if(cube26.is_cond_2_satisfied()) {
		  if(cube6.is_cond_4_satisfied()) ptl.push_back(ivec3(x,y,z));
		}
	      }
	    }
	  }
    }
  }

  int volume8::thin_apply(vector<ivec3> &ptl, thintype_t thintype) {
    int modified = 0;
    uint8 fg = 255;
    cube3x3 cube26; cube3x3_18 cube6;
    for(size_t i = 0; i < ptl.size(); i++) {
      int x = ptl[i][0], y = ptl[i][1], z = ptl[i][2];
      if(!is_finished(thintype, x,y,z,fg)) {
	// Fill cube for unfinished points. 
	for(int dx = -1; dx <= 1; dx++)
	  for(int dy = -1; dy <= 1; dy++)
	    for(int dz = -1; dz <= 1; dz++) {
	      cube26(dx,dy,dz) = v(x+dx,y+dy,z+dz);
	      cube6(dx,dy,dz) = v(x+dx,y+dy,z+dz);
	    }
	if(cube26.is_cond_2_satisfied()) {
	  if(cube6.is_cond_4_satisfied()) {
	    v(ptl[i]) = 0; modified++; 
	  }
	}
      }
    }
    return modified;
  }

  int volume8::thin_subiter(thindir_t direction, thintype_t thintype) {
    vector<ivec3> ptl; thinptl(ptl, direction, thintype);
    return thin_apply(ptl, thintype);
  }
  // Clear a 1 voxel border.
  void volume8::thin_init() {
    for(int x = 0; x < width; x++)
      for(int y = 0; y < height; y++) 
	for(int z = 0; z < depth; z++) {
	  if(x == 0 || y == 0 || z == 0) v(x,y,z) = 0;
	  if(x == width - 1 || y == height - 1 || z == depth - 1) v(x,y,z) = 0;
	}
  }

  void volume8::thin(thintype_t thintype) {
    thin_init();
    // Applying thinning subiterations.
    int modified = 0;
    do {
      modified = 0;
      // Apply for each of 6 directions.
      for(int d = 0; d < 6; d++) modified += thin_subiter((thindir_t)d, thintype);
    } while(modified > 0);
  }
  
  // Update histogram and cumulative value.
  void volume8::addhist(int hist[256], int med, int &cum,
			vector<ivec3> &s, int x, int y, int z) {
    for(size_t i = 0; i < s.size(); i++) {
      int val = v(x + s[i].x(), y + s[i].y(), z + s[i].z());
      hist[val]++;
      if(val < med) cum++;
    }
  }
  // Subtract out values from histogram and update cumulative value.
  void volume8::subhist(int hist[256], int med, int &cum,
			vector<ivec3> &s, int x, int y, int z) {
    for(size_t i = 0; i < s.size(); i++) {
      int val = v(x + s[i].x(), y + s[i].y(), z + s[i].z());
      hist[val]--;
      if(val < med) cum--;
    }
  }
  // Median algorithm that uses ideas from Huang et al
  void volume8::rankfilt(volume8 &dst, int r, float alpha, int zstart, int zstep) {
    //cout << "volume8:median zstart=" << zstart << " zstep=" << zstep << endl;
    vector<ivec3> s, sL, sR, sU, sD, sT, sB; 
    sphere_strel(s, r); // entire sphere centered at (0,0,0)
    // Prepare structural elements for updating fg value count.
    sphere_strel(sL, r,  1,  0, 0);     // remove sphere to right
    sphere_strel(sR, r, -1,  0, 0);     // remove sphere to left
    sphere_strel(sU, r,  0,  1, 0);     // remove sphere from below
    sphere_strel(sD, r,  0,  -1, 0);    // remove sphere from above
    sphere_strel(sT, r,  0, 0, zstep);  // remove sphere from top
    sphere_strel(sB, r,  0, 0, -zstep); // remove sphere from bottom

    int nalpha = max(int(s.size() * alpha), 1);  // NOTE: If nalpha = 0 then this gets messed up.
    int x = r, y = r, z = r + zstart, xstep = 1, ystep = 1; 
    int hist[256]; for(int i = 0; i < 256; i++) hist[i] = 0;

    int med = 0, cum = 0; addhist(hist, med, cum, s, x, y, z); // Fill initial histogram with entire sphere.

    med = cum = 0; while(med < 255 && cum < nalpha) cum += hist[med++];
    //cout << "med = " << med << " cum = " << cum << "/" << nalpha << endl;

    int d = depth - r, w = width - r, h = height - r;
    while(z < d) {
      while(r <= y && y < h) {
	while(r <= x && x < w) {
	  // Compute median.
	  if(cum < nalpha) { while(med < 255 && cum < nalpha) cum += hist[med++]; }
	  else { while(med > 0 && cum >= nalpha) { med--; cum -= hist[med]; } }
	  dst(x-r,y-r,z-r) = med;
	  x += xstep;
	  if((xstep > 0 && x < w) || (xstep < 0 && x >= r)) {
	    // Depending on direction of scan, pick correct STREL to
	    // use for updating fg counts.
	    vector<ivec3> &s1 = xstep > 0 ? sL : sR;
	    subhist(hist, med, cum, s1, x-xstep,y,z);
	    vector<ivec3> &s2 = xstep > 0 ? sR : sL;
	    addhist(hist, med, cum, s2, x,y,z);
	  }
	}
	xstep = -xstep; x += xstep; // flip direction of x-scan
	y += ystep;
	if((ystep > 0 && y < h) || (ystep < 0 && y >= r)) {
	  // Depending on the direction of the scan pick the corresponding
	  // strel for updating the y position.
	  vector<ivec3> &s1 = ystep > 0 ? sU : sD;
	  subhist(hist, med, cum, s1, x,y-ystep,z);
	  vector<ivec3> &s2 = ystep > 0 ? sD : sU;
	  addhist(hist, med, cum, s2, x,y,z);
	}
      }
      ystep = -ystep; y += ystep; // flip direction of y-scan.
      z += zstep;
      if( z < d ) {
	// The z-scan always goes deeper into the image volume.
	subhist(hist, med, cum, sT, x,y,z-zstep);
	addhist(hist, med, cum, sB, x,y,z);
      }
    }      
  }

  void volume8::fill_holes(uint8 fg, uint8 bg) {
    vector< vector<ivec3> > holes;  
    components(holes, 0, bg);
    if(holes.size() > 0) {
      size_t max_h = 0, max_vol = holes[0].size();
      for(size_t h = 1; h < holes.size(); h++) {
	if(holes[h].size() > max_vol) { max_vol = holes[h].size(); max_h = h; }
      }
      for(size_t h = 0; h < holes.size(); h++) {
	if(h == max_h) continue; // Skip largest.
	for(size_t i = 0; i < holes[h].size(); i++) v(holes[h][i]) = fg;
      }
    }
  }

  inline bool volume8::RemoveFT(int d2u, int d2v, int d2w, int ud, int vd, int wd) {
    int a = vd - ud, b = wd - vd, c = wd - ud;
    return c * d2v - b * d2u - a * d2w - a * b * c > 0;
  }
  void volume8::VoronoiFT(volumeT<ivec3> &F, int d, ivec3 j) {
    int maxnd = max(width, max(height, depth));
    vector<ivec3> g(maxnd+1);
    vector<int> h2(maxnd+1), g2(maxnd+1);
    int l = 0, nd = 0;
    if(d == 1) nd = width; else if(d == 2) nd = height; else if(d == 3) nd = depth;
    ivec3 xi = j;
    for(int i = 0; i < nd; i++) {
      xi[d-1] = i;
      ivec3 fi = F(xi); int dsq = distance3sq(xi, fi);
      if( fi[0] != -1 ) { // nil value
	if(l < 2) {  l = l+1; g[l] = fi; g2[l] = dsq; h2[l] = i; }
	else {
	  while( l >= 2 && RemoveFT(g2[l-1], g2[l], dsq, h2[l-1], h2[l], i) ) l = l - 1;
	  l = l + 1; g[l] = fi; g2[l] = dsq; h2[l] = i;
	}
      }
    }
    int ns = l; if(ns == 0) return;
    l = 1;
    for(int i = 0; i < nd; i++) {
      xi[d-1] = i;
      while( l < ns && g2[l] + sqr(h2[l] - i) > g2[l+1] + sqr(h2[l+1] - i) ) l = l + 1;
      F(xi) = g[l];
    }
  }
  void volume8::ComputeFT(volumeT<ivec3> &F, uint8 fg, int d, ivec3 j) {
    if(d == 1) {
      for(int x = 0; x < width; x++) {
	// fg values point to them selves
	if(v(x, j[1], j[2]) == fg) F(x, j[1], j[2]) = ivec3(x, j[1], j[2]);
	else F(x, j[1], j[2]) = ivec3(-1,-1,-1); // nil value
      }
    }
    else {
      if(d == 2) { for(int y = 0; y < height; y++) ComputeFT(F, fg, d - 1, ivec3(-1, y, j[2])); }
      else if(d == 3) { for(int z = 0; z < depth;  z++) ComputeFT(F, fg, d - 1, ivec3(-1,-1,z) ); }
    }
    if(d == 1) VoronoiFT(F, d, ivec3(-1, j[1], j[2]));
    else if(d == 2) {
      for(int x = 0; x < width; x++) VoronoiFT(F, d, ivec3(x, -1, j[2]));
    }
    else if(d == 3) {
      for(int x = 0; x < width; x++)
	for(int y = 0; y < height; y++) VoronoiFT(F, d, ivec3(x,y,-1));
    }
  }


  // From Vincent, L., "Morphological Grayscale Reconstruction in
  // Image Analysis: Applications and Efficient Algorithms," IEEE
  // Transactions on Image Processing, Vol. 2, No. 2, April, 1993,
  // pp. 176-201.  
  // J subset eq I: J can only have a fg value where I has a fg value!!!
  // J is the marker image
  // I is the input binary image.
  // output is in J
  void binary_reconstruct(volume8 &J, volume8 &I, uint8 fg, bool use_nhood26) {
    vector<ivec3> adj; adj.reserve(26);
    queue<ivec3> fifo;

    // Initialization of the queue with contour voxels of the marker image J
    for(int z = 0; z < I.depth; z++) 
      for(int y = 0; y < I.height; y++) 
	for(int x = 0; x < I.width; x++) {
	  ivec3 p(x,y,z);
	  if(J(p) == fg) {
	    bool add2fifo = false;
	    if(use_nhood26) I.nhood26(adj, p); else I.nhood6(adj, p);
	    for(size_t a = 0; a < adj.size(); a++) {
	      ivec3 q = adj[a];
	      if(J(q) != fg && I(p) == fg) { add2fifo = true; break; }
	    }
	    if(add2fifo) fifo.push(p);
	  }
	}
    
    // Propagation.
    while(!fifo.empty()) {
      ivec3 p = fifo.front(); fifo.pop();
      if(use_nhood26) I.nhood26(adj, p); else I.nhood6(adj, p);
      for(size_t a = 0; a < adj.size(); a++) {
	ivec3 q = adj[a];
	if(J(q) != fg && I(q) == fg) {
	  J(q) = fg;
	  fifo.push(q);
	}
      }
    }
  }

  void check_IJ(volume8 &I, volume8 &J) {
    for(int z = 0; z < I.depth; z++) 
      for(int y = 0; y < I.height; y++) 
	for(int x = 0; x < I.width; x++) {
	  if(! (J(x,y,z) <= I(x,y,z) ) ) {
	    cerr << "check_IJ failed" << endl; 
	    exit(1);
	  }
	}    
  }

  void reconstruct_raster(volume8 &J, volume8 &I, bool use_nhood26) {
    vector<ivec3> adj, check; adj.reserve(26); check.reserve(26);

    check_IJ(I,J); // confirm J <= I (assumes fg = 255).

    int updated;
    do {
      updated = 0;

      // Raster order.
      for(int z = 0; z < I.depth; z++)
	for(int y = 0; y < I.height; y++) 
	  for(int x = 0; x < I.width; x++) {
	    ivec3 p(x,y,z);
	    if(use_nhood26) I.nhood26(adj, p); else I.nhood6(adj, p);
	    check.clear();
	    // Get N+ set, neighbors reached before p during raster scanning (left to right top to bottom)
	    for(size_t a = 0; a < adj.size(); a++) {
	      if(adj[a].x() < x || adj[a].y() < y || adj[a].z() < z) check.push_back(adj[a]);
	    }
	    int maxJ = J(p);
	    for(size_t c = 0; c < check.size(); c++) maxJ = max((int)J(check[c]), maxJ);
	    int newJ = min(maxJ, (int)I(p));
	    if(J(p) != newJ) { J(p) = newJ; updated++; }
	  }

      // Anti-raster order
      for(int z = I.depth - 1; z >= 0; z--) 
	for(int y = I.height - 1; y >= 0; y--) 
	  for(int x = I.width - 1; x >= 0; x--) {
	    ivec3 p(x,y,z);
	    if(use_nhood26) I.nhood26(adj, p); else I.nhood6(adj, p);
	    check.clear();
	    // Get N- set, neighbors reached *after* p during raster scanning (left to right top to bottom)
	    for(size_t a = 0; a < adj.size(); a++) {
	      if(adj[a].x() > x || adj[a].y() > y || adj[a].z() > z) check.push_back(adj[a]);
	    }
	    int maxJ = J(p);
	    for(size_t c = 0; c < check.size(); c++) maxJ = max((int)J(check[c]), maxJ);
	    int newJ = min(maxJ, (int)I(p));
	    if(J(p) != newJ) { J(p) = newJ; updated++; }
	  }

    } while( updated > 0) ;

  }

  
  // Hybrid Algorithm from Vincent, L., "Morphological Grayscale
  // Reconstruction in Image Analysis: Applications and Efficient
  // Algorithms," IEEE Transactions on Image Processing, Vol. 2,
  // No. 2, April, 1993, pp. 176-201.
  void grayscale_reconstruct(volume32 &J, volume32 &I, bool use_nhood26) {
    vector<ivec3> adj, check; adj.reserve(26); check.reserve(26);
    queue<ivec3> fifo;
    for(int z = 0; z < I.depth; z++) 
      for(int y = 0; y < I.height; y++) 
	for(int x = 0; x < I.width; x++) {
	  if(! (J(x,y,z) <= I(x,y,z) ) ) {
	    cerr << "grayscale_reconstruct(): check_IJ failed" << endl; 
	    exit(1);
	  }
	}    
    // Scan Di in raster order.
    for(int z = 0; z < I.depth; z++)
      for(int y = 0; y < I.height; y++) 
	for(int x = 0; x < I.width; x++) {
	    ivec3 p(x,y,z);
	    if(use_nhood26) I.nhood26(adj, p); else I.nhood6(adj, p);
	    check.clear();
	    // Get N+ set, neighbors reached before p during raster scanning (left to right top to bottom)
	    for(size_t a = 0; a < adj.size(); a++) {
	      if(adj[a].x() < x || adj[a].y() < y || adj[a].z() < z) check.push_back(adj[a]);
	    }
	    int maxJ = J(p);
	    for(size_t c = 0; c < check.size(); c++) maxJ = max((int)J(check[c]), maxJ);
	    int newJ = min(maxJ, (int)I(p));
	    J(p) = newJ; 
	  }
    
    // Anti-raster order
    for(int z = I.depth - 1; z >= 0; z--) 
      for(int y = I.height - 1; y >= 0; y--) 
	for(int x = I.width - 1; x >= 0; x--) {
	  ivec3 p(x,y,z);
	  if(use_nhood26) I.nhood26(adj, p); else I.nhood6(adj, p);
	  check.clear();
	  // Get N- set, neighbors reached *after* p during raster scanning (left to right top to bottom)
	  for(size_t a = 0; a < adj.size(); a++) {
	    if(adj[a].x() > x || adj[a].y() > y || adj[a].z() > z) check.push_back(adj[a]);
	  }
	  int maxJ = J(p);
	  for(size_t c = 0; c < check.size(); c++) maxJ = max((int)J(check[c]), maxJ);
	  int newJ = min(maxJ, (int)I(p));
	  J(p) = newJ;
	  for(size_t c = 0; c < check.size(); c++) {
	    ivec3 q = check[c];
	    if(J(q) < J(p) && J(q) < I(q)) fifo.push(p); 
	  }
	}

    while(!fifo.empty()) {
      ivec3 p = fifo.front(); fifo.pop();
      if(use_nhood26) I.nhood26(adj, p); else I.nhood6(adj, p);
      for(size_t a = 0; a < adj.size(); a++) {
	ivec3 q = adj[a];
	if( J(q) < J(p) && I(q) != J(q) ) {
	  J(q) = min(J(q), I(q));
	  fifo.push(q);
	}
      }
    }

  }

  // Computes shortest path between two points in a volume.
  float ShortestPath(vector<ivec3> &path, volume8 &v, ivec3 src, ivec3 targ, uint8 fg) {
    path.clear();
    // Check source and target positions 
    if(v(src) != fg || v(targ) != fg) return 0;  

    // Create mapping between index <-> position.
    vector<ivec3> idx2pos;
    volume32 pos2idx(v.width, v.height, v.depth); pos2idx.fill(-1); 
    int idx = 0;
    for(int z = 0; z < v.depth; z++) 
      for(int y = 0; y < v.height; y++) 
	for(int x = 0; x < v.width; x++) {
	  if(v(x,y,z) == fg) {
	    ivec3 pos(x,y,z);
	    idx2pos.push_back(pos);
	    pos2idx(pos) = idx;
	    idx++;
	  }
	}
    
    int N = idx2pos.size();
    const float inf = numeric_limits<float>().max();

    vector<float> dist(N);
    util::PQi<float> Q(N, dist);
    vector<int> previous(N);
    for(idx = 0; idx < N; idx++) { dist[idx] = inf;  previous[idx] = -1; }

    dist[pos2idx(src)] = 0;     // Fill priority queue.
    for(int idx = 0; idx < N; idx++) Q.insert(idx);

    vector<ivec3> adjtest;
    while(!Q.empty()) {
      int u = Q.extractMin();
      if(dist.at(u) == inf) break;  //?

      ivec3 &pos = idx2pos[u];
      v.nhood26(adjtest, pos);
      for(size_t a = 0; a < adjtest.size(); a++) {
	int v = pos2idx(adjtest[a]);
	if(v == -1) continue; // Skip bg adjacencies.
	float alt = dist[u] + distance3sq(idx2pos[v],idx2pos[u]);
	if(alt < dist[v]) {
	  dist[v] = alt; Q.decreaseKey(v);
	  previous[v] = u;
	}
      }
    }
    int u = pos2idx(targ);
    if(previous[u] < 0) return 0; // no path

    while(previous[u] >= 0) {
      path.push_back(idx2pos[u]);
      u = previous[u];
    }
    path.push_back(idx2pos[u]);
    reverse(path.begin(), path.end());
    return dist[pos2idx(targ)];
  }

  // Generates a synthetic volume using a feature transform.
  void GenSynth(volume8 &v, int num_cells, float /* max_dist */) {
    volume32 labels(v.width, v.height, v.depth);
    v.fill(0);
    labels.fill(0);
    srand48(time(NULL));

    for(int label = 1; label <= num_cells; label++) {
      // Sample two points in volume and create a vector between them.
      int x = drand48() * v.width, y = drand48() * v.height, z = drand48() * v.depth;
      v(x,y,z) = 255;
      labels(x,y,z) = label;
    }

    // Creat an fg volume and a corresponding label volume.
    /*for(int label = 1; label <= num_cells; label++) {
      // Sample two points in volume and create a vector between them.
      int x1 = drand48() * v.width, y1 = drand48() * v.height, z1 = drand48() * v.depth;
      int x2 = drand48() * v.width, y2 = drand48() * v.height, z2 = drand48() * v.depth;
      vec3 p1(x1, y1, z1), p2(x2, y2, z2);
      vec3 dir = p2 - p1;
      dir.normalize();

      float dthresh = max_dist * (1.0 - drand48());

      vector<ivec3> isect; 
      v.ray_intersect(isect, p1, dir, dthresh);
      if(isect.size() == 0) cout << "isect.size() == 0" << endl;
      for(size_t r = 0; r < isect.size(); r++) {
	v(isect[r]) = 255;
	labels(isect[r]) = label;
      }
      }*/
    // 
    volumeT<ivec3> F(v.width, v.height, v.depth);
    cout << "computing FT" << endl;
    v.ComputeFT(F, 255);
    cout << "DONE computing FT" << endl;
    for(int z = 0; z < F.depth; z++) 
      for(int y = 0; y < F.height; y++) 
	for(int x = 0; x < F.width; x++) {
	  // Use FT to get nearest labeled voxel.
	  int label = labels(F(x,y,z));
	  if(label == 0) cout << "FT FAILED ********" << endl;
	  // Assign curent voxel the label.
	  labels(x,y,z) = label;
	}
    v.fill(0);
    for(int z = 0; z < v.depth; z++) 
      for(int y = 0; y < v.height; y++) 
	for(int x = 0; x < F.width; x++) {
	  // Check for a boundary.
	  for(int dz = -1; dz <= 1; dz++)
	    for(int dy = -1; dy <= 1; dy++)
	      for(int dx = -1; dx <= 1; dx++) {
		if(dx == 0 && dy == 0 && dz == 0) continue;
		int x2 = x+dx, y2 = y+dy, z2 = z+dz;
		if(x2 < 0 || y2 < 0 || z2 < 0) continue;
		if(x2 >= v.width || y2 >= v.height || z2 >= v.depth) continue;

		if(labels(x2,y2,z2) != labels(x,y,z)) {
		  v(x,y,z) = 255;
		}
	      }
	}
  }

  // Counts internal regions. You would find a lot of these if two
  // cells with thinned boundaries were merged. Works quite well,
  // leads to some oversegmentation.
  int HS2_CheckComp(vector<ivec3> &comp) {
    volume8 v(comp, 2);

    // Dilate region.
    vector<ivec3> dvox; v.dilate(dvox, 1, 255, 0);
    for(size_t i = 0; i < dvox.size(); i++) v(dvox[i]) = 255;

    // Erode.. gets rid of "small" inward regions on surface 
    // addresses noise from thinning.
    vector<ivec3> erd;  v.dilate(erd, 1, 0, 255);
    for(size_t i = 0; i < erd.size(); i++) v(erd[i]) = 0;

    // Find dilated voxels that are adjacent to the "outside" of
    // current region.
    vector<ivec3> adj(26);
    int internal_cnt = 0;
    for(size_t i = 0; i < dvox.size(); i++) {
      bool dvox_outside = false;
      // outside if bg value after eroding
      if(v(dvox[i]) == 0) dvox_outside = true;
      else { // or a neighbor is a bg value
	v.nhood26(adj, dvox[i]);
	for(size_t a = 0; a < adj.size(); a++) {
	  if(v(adj[a]) == 0) { dvox_outside = true; break; }
	}
      }
      // If dilated region isn't a bg value or neighboring, it is internal.
      if(!dvox_outside) internal_cnt++;
    }

    return internal_cnt;
  }

  // Clear out internal voxels from a region.
  void HS2_Clear(volume8 &orig, vector<ivec3> &comp) {
    volume8 v(comp, 2);

    // Dilate region.
    vector<ivec3> dvox; v.dilate(dvox, 1, 255, 0);
    for(size_t i = 0; i < dvox.size(); i++) v(dvox[i]) = 255;

    // Erode.. gets rid of "small" inward regions on surface 
    // addresses noise from thinning.
    vector<ivec3> erd;  v.dilate(erd, 1, 0, 255);
    for(size_t i = 0; i < erd.size(); i++) v(erd[i]) = 0;

    // Find dilated voxels that are adjacent to the "outside" of
    // current region.
    vector<ivec3> adj(26);
    //int internal_cnt = 0;
    for(size_t i = 0; i < dvox.size(); i++) {
      bool dvox_outside = false;
      // outside if bg value after eroding
      if(v(dvox[i]) == 0) dvox_outside = true;
      else { // or a neighbor is a bg value
	v.nhood26(adj, dvox[i]);
	for(size_t a = 0; a < adj.size(); a++) {
	  if(v(adj[a]) == 0) { dvox_outside = true; break; }
	}
      }
      // If dilated region isn't a bg value or neighboring, it is internal.
      if(!dvox_outside) {
	orig(dvox[i][0] + v.x0, dvox[i][1] + v.y0, dvox[i][2] + v.z0) = 0;
	ivec3 missing(dvox[i][0] + v.x0, dvox[i][1] + v.y0, dvox[i][2] + v.z0);
	comp.push_back(missing);
      }
    }

  }

  void HS2_Clear(vector<ivec3> &comp) {
    volume8 v(comp, 2);

    // Dilate region.
    vector<ivec3> dvox; v.dilate(dvox, 1, 255, 0);
    for(size_t i = 0; i < dvox.size(); i++) v(dvox[i]) = 255;

    // Erode.. gets rid of "small" inward regions on surface 
    // addresses noise from thinning.
    vector<ivec3> erd;  v.dilate(erd, 1, 0, 255);
    for(size_t i = 0; i < erd.size(); i++) v(erd[i]) = 0;

    // Find dilated voxels that are adjacent to the "outside" of
    // current region.
    vector<ivec3> adj(26);
    //int internal_cnt = 0;
    for(size_t i = 0; i < dvox.size(); i++) {
      bool dvox_outside = false;
      // outside if bg value after eroding
      if(v(dvox[i]) == 0) dvox_outside = true;
      else { // or a neighbor is a bg value
	v.nhood26(adj, dvox[i]);
	for(size_t a = 0; a < adj.size(); a++) {
	  if(v(adj[a]) == 0) { dvox_outside = true; break; }
	}
      }
      // If dilated region isn't a bg value or neighboring, it is internal.
      if(!dvox_outside) {
	ivec3 missing(dvox[i][0] + v.x0, dvox[i][1] + v.y0, dvox[i][2] + v.z0);
	comp.push_back(missing);
      }
    }

  }

  bool EDT_Check(vector<ivec3> &merge_comp, vector<ivec3> &c1, vector<ivec3> &c2) {
    volume8 merged(merge_comp);
    vector<ivec3> erd; 
    merged.dilate(erd, 2, 0, 255);

    // Erode to break components.
    for(size_t e = 0; e < erd.size(); e++) merged(erd[e]) = 0;
    vector< vector<ivec3> > comps;
    merged.components(comps, 0, 255);

    // If two or more components...
    if(comps.size() >= 2) {
      volume32 labels(merged.width, merged.height, merged.depth);
      ivec3 p(merged.x0, merged.y0, merged.z0);

      for(size_t i = 0; i < c1.size(); i++) labels(c1[i] - p) = 1;
      for(size_t i = 0; i < c2.size(); i++) labels(c2[i] - p) = 2;

      // Check to see if broken components span two regions.
      // If component does, then we're OK to merge.
      for(size_t c = 0; c < comps.size(); c++) {
	int label_idx = -1;
	vector<ivec3> &comp = comps[c];
	for(size_t i = 0; i < comp.size(); i++) {
	  if(labels(comp[i]) == 0) continue;
	  if(label_idx < 0) label_idx = labels(comp[i]);
	  else if(label_idx != labels(comp[i])) return true;
	}
      }
      return false;
    }
    else return true;
  }


  void SeedSegment(volume8 &edge_vol,
		   vector< vector<ivec3> > &res, volume8 &v, vector< vector<ivec3> > &seeds, int min_comp, int max_comp) {
    if(seeds.size() == 0) return; // Nothing to seed, quit.
    int w = v.width, h = v.height, d = v.depth;
    volume32 labels(w,h,d), negEDT(w,h,d);
    
    negEDT.ComputeEDT(v);  // Compute negated EDT. TODO: Move this out and/or use multithreading.
    negEDT.negate();

    for(int z = 0; z < d; z++)
      for(int y = 0; y < h; y++) 
	for(int x = 0; x < w; x++) negEDT(x,y,z) -= 255 - (int)edge_vol(x,y,z);


    labels.seed(seeds);  // Initialized seeds.
    labels.seed(v, 255, -1);  // Exclude edge voxels

    labels.watershed(negEDT);   // Propagate labels using given negated EDT.
    vector< vector<ivec3> > cur; 
    // Find components with a label > 0, keeps empty components.
    labels.components(cur, false);  

    // Keep those in size range.
    for(size_t c = 0; c < cur.size(); c++) {
      if(min_comp <= (int)cur[c].size() && (int)cur[c].size() <= max_comp) res.push_back(cur[c]);
    }
  }



  bool HS2_BadCentroid(vector<ivec3> &comp) {
    vec3 cent_v = geom::Centroid3(comp);
    ivec3 cent(cent_v[0], cent_v[1], cent_v[2]);

    
    // Make sure centroid is inside the cell.
    volume8 filled(comp, 2);

    cent[0] -= filled.x0;
    cent[1] -= filled.y0;
    cent[2] -= filled.z0;

    if(filled(cent) == 0)  return true;

    return false;
  }

  struct adj_t {  
    int idx;
    float contact_score;
    adj_t() { idx = -1; contact_score = -1; }
    adj_t(int idx, float contact_score) : idx(idx), contact_score(contact_score) { }
    bool operator < (const adj_t &b) const { return contact_score > b.contact_score; }
  };

  // NOTE: HMerge() primarily seems to address the effects of salt-pepper noise in 
  // edge volume. 
  void HMerge_Nucs(volume32 &labels, volume32 *nucs, int internal_cnt,  
		   vector< vector<ivec3> > &final, int max_comp, vector< vector<ivec3> > &nuc_vox) {
    int w = labels.width, h = labels.height, d = labels.depth;
    vector<bool> failed_merge(final.size(), false);
    labels.seed(final);

    int max_sz = 0;
    for(int fidx = 0; fidx < (int)final.size(); fidx++) { if((int)final[fidx].size() > max_sz) max_sz = final[fidx].size(); }

    int merge_cnt = 0;
    while(true) {
      // Get smallest region w/o merge failure or hasn't been merged already.
      int min_fidx = -1, min_fsz = w * h * d;
      for(int fidx = 0; fidx < (int)final.size(); fidx++) {
	if(final[fidx].size() == 0) continue; // Region has been merged.
	if(failed_merge[fidx]) continue;      // No more merges possible w/o failure.
	if((int)final[fidx].size() < min_fsz) {
	  min_fsz = final[fidx].size();
	  min_fidx = fidx;
	}
      }
      if(min_fidx == -1) break; // Nothing found, stop.
      if(min_fsz > max_comp) break; // Stop after smallest component reaches min_comp size.

      // Get adjacencies and their size.
      vector<ivec3> &smallest = final[min_fidx], hood(6);
      vector<adj_t> adj;
      for(size_t s = 0; s < smallest.size(); s++) {
	// Look in 6-hood for an adjacency.
	labels.nhood6(hood, smallest[s]);
	for(size_t h = 0; h < hood.size(); h++) {
	  int adj_idx = labels(hood[h]) - 1;
	  if(adj_idx < 0 || adj_idx == min_fidx) continue; // Don't self merge. 
	  if(failed_merge.at(adj_idx)) continue; // Skip merge to a failed merge.
	  // Look for adjacency and increment contact area.
	  bool found_adj = false;
	  for(size_t a = 0; a < adj.size(); a++) {
	    if(adj[a].idx == adj_idx) {
	      // Smallest contact area, inversely weighted by volume,
	      // favor merging to smaller neighbors.
	      adj[a].contact_score += 1.0f / float(final[adj_idx].size());
	      found_adj = true;
	      break;
	    }
	  }
	  if(!found_adj) adj.push_back(adj_t(adj_idx, 1.0f / float(final[adj_idx].size()) ));
	}
      }
      // No adjacencies, stop trying to merge.
      if(adj.size() == 0) { failed_merge[min_fidx] = true; continue; }

      // Get region with largest weighted contact area. 
      sort(adj.begin(), adj.end()); 
      int a = 0;
      int adj_idx = adj[a].idx;
      if(failed_merge[adj_idx]) { cerr << "failed merge shouldn't happen" << endl; exit(1);  }
      if(final[adj_idx].size() == 0) { cerr << "shouldn't happen" << endl; exit(1); }

      vector<ivec3> test_merge = final.at(adj_idx);
      for(size_t i = 0; i < smallest.size(); i++) test_merge.push_back(smallest[i]);
      
      if((int)test_merge.size() > max_comp) { failed_merge[adj_idx] = true; continue; }

      // TODO: Check for multiple nucs.
      if(nucs != NULL) {
	vector<int> nuc_idxs, overlap_size;
	for(size_t i = 0; i < test_merge.size(); i++) {
	  int nuc_idx = nucs->v(test_merge[i]) - 1;
	  if(nuc_idx < 0) continue;
	  bool found_idx = false;
	  for(size_t k = 0; k < nuc_idxs.size(); k++) { 
	    if(nuc_idx == nuc_idxs[k]) { 
	      found_idx = true; 
	      overlap_size.at(k)++;
	      break; 
	    } 
	  }
	  if(!found_idx) {
	    nuc_idxs.push_back(nuc_idx);
	    overlap_size.push_back(1);
	  }
	}
	if(nuc_idxs.size() > 1) { 
	  int failed_cnt = 0;
	  for(size_t k = 0; k < nuc_idxs.size(); k++) {
	    if(float(overlap_size.at(k)) / float( nuc_vox.at( nuc_idxs.at(k) ).size() ) > 0.1) failed_cnt++; 
	  }
	  if(failed_cnt > 1) {
	    failed_merge[adj_idx] = true;
	    continue; 
	  }
	}
      }

      int HS2_check = HS2_CheckComp(test_merge);
      //cout << "HS2_check=" << HS2_check << endl;
      if(HS2_check > internal_cnt) { failed_merge[adj_idx] = true; continue; }

      if(HS2_BadCentroid(test_merge)) { failed_merge[adj_idx] = true; continue; }

      if((int)test_merge.size() >  max_sz / 2) {
	if(EDT_Check(test_merge, smallest, final.at(adj_idx)) == false) {
	  failed_merge[adj_idx] = true; 
	  continue; 
	}
      }

      // Actually perform the merge.
      //if(test_merge.size() < max_sz / 4) 
      HS2_Clear(test_merge); // Clears out internal voxels.

      for(size_t i = 0; i < smallest.size(); i++)  // Re-label previous small region.
	labels(smallest[i]) = adj_idx+1;

      // Update region.
      final.at(adj_idx) = test_merge;	  

      merge_cnt++;
      //cout << "smallest.size()=" << smallest.size() << endl;
      //merge_sz.at(adj_idx) += merge_sz.at(adj_idx) + (int)smallest.size();

      smallest.clear(); // Clear old region.
    }
    // Debugging to visualize merged cells.
    //for(size_t m = 0; m < merge_sz.size(); m++) { if(merge_sz[m] <= 1000) final[m].clear();  }

    cout << "merge_cnt=" << merge_cnt << endl;
  }
  

  // Approach has a tendency to oversegment. 
  void HSegment3_Nucs(volume8 &edge_vol, volume8 &vorig, 
  		      vector< vector<ivec3> > &nuc_vox, 
  		      int dmin, int dmax, 
  		      int internal_cnt,  
  		      vector< vector<ivec3> > &final, 
  		      int min_comp, int max_comp, int noise_comp) {
    int w = vorig.width, h = vorig.height, d = vorig.depth;
    vector< vector<ivec3> > init_seeds, seeds_cleaned, cur, checked;
    volume8 v(w,h,d);

    v.copy(vorig); // Keep a copy of the original volume.

    // Use Nuclei to seed Watershed segmentation.
    if(nuc_vox.size() > 0)  {
      cout << "nuc_vox.size()=" << nuc_vox.size() << endl;
      cur.clear(); 
      SeedSegment(edge_vol, cur, v, nuc_vox, min_comp, max_comp);
      for(size_t c = 0; c < cur.size(); c++) {
	vector<ivec3> &comp = cur[c];
	// Check internal voxel count.
	if(HS2_CheckComp(comp) > internal_cnt) continue; // Too many internal voxels.

	final.push_back(comp); // OK, keep as a final component.
      }
    }
    cout << "[seed] final.size()=" << final.size() << endl;

    volume32 labels(w,h,d), negEDT(w,h,d), negEDT_water(w,h,d);
    volume32 *nucs = NULL;
    if(nuc_vox.size() > 0) nucs = new volume32(w,h,d);

    volume8 markers(w,h,d);

    // Eliminate finalized regions from watershed.
    for(size_t f = 0; f < final.size(); f++) 
      for(size_t i = 0; i < final[f].size(); i++) v(final[f][i]) = 255;   

    if(nuc_vox.size() > 0) nucs->seed(nuc_vox); // Seed nuclei regions to allow for checking.

    negEDT.ComputeEDT(v); negEDT.negate(); // Compute negated EDT. 

    for(int dthresh = dmin; dthresh <= dmax; dthresh++) {
      cout << "dthresh=" << dthresh << endl;
      negEDT.below_threshold(markers, -(dthresh * dthresh));

      // Find seed regions.
      init_seeds.clear(); 
      markers.components(init_seeds, 0, 255);
      if(init_seeds.size() == 0) break; // No seeds, quit.

      // NOTE: This seed cleaning step limits the size of components found to noise_comp voxels.
      labels.seed(init_seeds);  // Get initial seed regions.
      labels.seed(v, 255, -1);  // Exclude from Watershed transform from edge and previously segmented voxels.

      negEDT_water.copy(negEDT);
      for(int z = 0; z < d; z++)
	for(int y = 0; y < h; y++) 
	  for(int x = 0; x < w; x++) negEDT_water(x,y,z) -= 255 - (int)edge_vol(x,y,z);
	    
      labels.watershed(negEDT_water); // Propagate labels using given negated EDT.

      cur.clear(); labels.components(cur, false);  // Find components with a label > 0.

      // Re-seed killing noise junk regions (smaller than noise_comp size).
      seeds_cleaned.clear();
      for(size_t c = 0; c < cur.size(); c++) {
  	if((int)cur[c].size() >= noise_comp) seeds_cleaned.push_back(cur[c]); 
      }

      labels.seed(seeds_cleaned); // Initialize new cleaned seed regions.
      labels.seed(v, 255, -1);    // Exclude from Watershed transform from edge and previously segmented voxels.
      labels.watershed(negEDT_water);   // Propagate labels using given negated EDT.

      cur.clear(); 
      labels.components(cur, false);  // Find components with a label > 0.

      for(size_t c = 0; c < cur.size(); c++) {
	vector<ivec3> &comp = cur[c];
	if((int)comp.size() > max_comp) continue;  // Region is too big (needs more splitting), skip.

	// Check internal voxel count.
	if(HS2_CheckComp(comp) > internal_cnt) continue; // Too many internal voxels, break some more.

	if(HS2_BadCentroid(comp)) continue;

	// Region doesn't require further breaking by EDT threshold, save as checked.
	checked.push_back(comp);
      }

      // Block out checked regions from further splitting.
      for(size_t f = 0; f < checked.size(); f++) {
	vector<ivec3> &comp = checked[f];
	for(size_t i = 0; i < comp.size(); i++) {
	  negEDT(comp[i]) = 0;  // Kill finished region in EDT
	  v(comp[i]) = 255;     // Eliminate it from Watershed regions for increasing dthresh.
	}
      }
      cout << "checked.size()=" << checked.size() << endl;
      // increase d-max and break regions further
    }

    for(size_t c = 0; c < checked.size(); c++)	final.push_back(checked[c]);

    // Get anything that's left.
    vector< vector<ivec3> > remaining;
    v.components(remaining, noise_comp, 0);
    cout << "remaining.size()=" << remaining.size() << endl;
    for(size_t r = 0; r < remaining.size(); r++) 
      if((int)remaining[r].size() <= max_comp) final.push_back(remaining[r]);


    sort(final.begin(), final.end(), cmp_size_gt<ivec3>()); 
    for(size_t f = 0; f < final.size(); f++) {
      HS2_Clear(vorig, final[f]);
    }

    // NOTE: use smaller radius for rank filtering 0.95 and 0.05, but using bigger sigma small and threhsold
    for(int iter = 0; iter < 2; iter++) {
      cout << "merge_iter=" << iter + 1 << endl;
      HMerge_Nucs(labels, nucs, internal_cnt, final, max_comp, nuc_vox);
    }

    int comp_cnt = 0;
    for(size_t f = 0; f < final.size(); f++)  if((int)final[f].size() >= min_comp) comp_cnt++;
    vector< vector<ivec3> > final2(comp_cnt);
    int f2pos = 0;
    for(size_t f = 0; f < final.size(); f++)  if((int)final[f].size() >= min_comp) final2[f2pos++] = final[f];
    final = final2;

    cout << "final.size()=" << final.size() << endl;
    sort(final.begin(), final.end(), cmp_size_gt<ivec3>()); 

    if(nucs != NULL) delete nucs;
  }


  ivec3 hausdorff_align(vector<ivec3> &A, vector<ivec3> &B) {
    volume8 A_ovol(A,2), B_ovol(B,2);
    vector<ivec3> A_outline, B_outline, AB_outline; 

    // Get voxel outlines.
    ivec3 dA_ovol(A_ovol.x0, A_ovol.y0, A_ovol.z0);
    ivec3 dB_ovol(B_ovol.x0, B_ovol.y0, B_ovol.z0);

    A_ovol.outline_vox(A_outline); B_ovol.outline_vox(B_outline);

    // Shift outlines back to non-volume frame.
    for(size_t i = 0; i < A_outline.size(); i++) {
      A_outline[i] += dA_ovol;
      AB_outline.push_back(A_outline[i]);
    }
    for(size_t i = 0; i < B_outline.size(); i++) {
      B_outline[i] += dB_ovol;
      AB_outline.push_back(B_outline[i]);
    }

    // Get combined volume.
    volume8 AB_ovol(AB_outline, 5);

    ivec3 dAB(AB_ovol.x0, AB_ovol.y0, AB_ovol.z0);
    for(size_t i = 0; i < A_outline.size(); i++) A_outline[i] -= dAB;
    for(size_t i = 0; i < B_outline.size(); i++) B_outline[i] -= dAB;

    AB_ovol.fill(0);
    AB_ovol.set(B_outline, 255);
    volume32 EDT(AB_ovol.width, AB_ovol.height, AB_ovol.depth);
    EDT.ComputeEDT(AB_ovol);

    ivec3 min_shift;
    EDT.frac_hausdorff_match(min_shift, A_outline, 5, 0.95);
    return min_shift;
  }


}; // namespace vol




