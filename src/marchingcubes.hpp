#ifndef __MARCHINGCUBES_HPP__

#include <vector>

#include "volume.hpp"
#include "geometry.hpp"

namespace marching {
  using namespace geom;
  using namespace vol;
  using namespace std;

  // c = cxcxc cube size Applies the Marching Cubes algorithm. Returns
  // output as a vertex + face table.  This output allows face/face
  // adjacencies and vertex normals to be computed.  See mesh.[h,c]pp.
  void MarchingCubes(float isolevel, volume8 &v, vector<vec3> &vtable, vector<face> &ftable, int c);
}

#endif // __MARCHINGCUBES_HPP__
