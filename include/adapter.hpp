#ifndef ADAPTER_HPP_
#define ADAPTER_HPP_

#include <string>

#include "model/base.hpp"

namespace FCA
{
class Adapter
{
public:
  Adapter(
    const FCA::FisheyeCameraModel * const input_model,
    FCA::FisheyeCameraModel * const output_model);

  void adapt();

  void compare_image(const std::string & image_path, bool show = false);

  void evaluate();
  
private:
  std::vector<Eigen::Vector2d> sample_points(int32_t width,int32_t height, int n);

  const FCA::FisheyeCameraModel * const input_model_;
  FCA::FisheyeCameraModel * const output_model_;
  std::vector<Eigen::Vector2d> point2d_vec_;
  std::vector<Eigen::Vector3d> point3d_vec_;
};
}  // namespace FCA
#endif
