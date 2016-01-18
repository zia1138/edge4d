#ifndef __PROCTHREADS__HPP__
#define __PROCTHREADS__HPP__


#include "volume.hpp"
#include "mesh.hpp"

enum gauss_dim { DimX = 0,  DimY = 1, DimZ = 2 };

using namespace vol;
// NOTE: On Mac OSX, cout/cin are not thread safe. They will cause
// crashes if you output in multithreaded code.

// Multithreaded processing steps.

// Multithreaded Gaussian blur.
void blur(volume8 *v, gauss_dim dim, float sigma, int num_threads = 1);
// Stretches intensity values to span [0,65K] range.
void stretch16(volume16 *v, int num_threads);
// Multithreaded histogram mapping. Adjustment or match to global histogrma.
void happly(volume8 *v, int nx, int ny, int nz, float fmin, float fmax, int num_threads = 1);
// Dilte or erode an image. 
void dilate(volume8 &v, int r, float alpha, uint8 bg, uint8 fg, bool invert, uint8 uval, int num_threads);
// Multithreaded adapt of meshes to fit volume intensity.
void refine(volume8 *v, float dist, int max_iter, float alpha, float stepsize,
	    vector<geom::mesh *> &cells, int num_threads = 1);
// Multithreaded rank filtering.
void rankfilt(volume8 &v, int r, float alpha, int num_threads);
// Multithreaded Euclidian distance transform. 
void compute_edt(volume32 &D, volume8 &I, uint8 fg, int num_threads);
// Multithreaded thinning. Thin a binary volume where fg = 255.
void thin(volume8 &v, thintype_t thintype, int num_threads);

// Multithreaded image scaling.  All of these either use trilinear
// interpolation (todo, use tricubic?) to scale down or blur prior to
// nearest nighbor scaling to avoid antialiasing.
volume8 *scaleToZ(volume8 *src,  float sxy,  int num_threads); // Scales (x,y) axes.
volume8 *scaleToXY(volume8 *src, float sz,   int num_threads); // Scales z axis.
volume8 *scaleUpXY(volume8 *src, float sx, float sy, int num_threads);
volume8 *scaleXYZ(volume8 *src,  float sxyz, int num_threads); // Scales (x,y,z) axes.


#endif // __PROCTHREADS__HPP__
