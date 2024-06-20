#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <string>

#include "model/base.hpp"

namespace FCA
{
template <typename T>
struct PointCloud
{
  struct Point
  {
    T x, y, z;
  };

  using coord_t = T;

  std::vector<Point> pts;

  inline size_t kdtree_get_point_count() const { return pts.size(); }

  inline T kdtree_get_pt(const size_t idx, const size_t dim) const
  {
    if (dim == 0)
      return pts[idx].x;
    else if (dim == 1)
      return pts[idx].y;
    else
      return pts[idx].z;
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX & /* bb */) const
  {
    return false;
  }
};

FisheyeCameraModelPtr Create(const std::string & model_name, const std::string & dataset_path);

}  // namespace FCA
#endif
