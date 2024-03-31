#ifndef ADAPTER_HPP_
#define ADAPTER_HPP_

#include <opencv2/opencv.hpp>
#include <opencv2/viz.hpp>
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

  void evaluate();

  void display_point2d_vec(
    const std::string & name, const std::vector<Eigen::Vector2d> & point2d_vec) const;

  void display_point3d_vec(
    const std::string & name, const std::vector<Eigen::Vector3d> & point3d_vec) const;

  cv::Mat recover_image(
    double width, double height, const std::vector<Eigen::Vector3d> & point3d_vec) const;

  void compare_image(const std::string & image_path, bool show = false);

private:
  std::vector<Eigen::Vector2d> sample_points(int32_t width, int32_t height, int n);

  const FCA::FisheyeCameraModel * const input_model_;
  FCA::FisheyeCameraModel * const output_model_;
  std::vector<Eigen::Vector2d> point2d_vec_;
  std::vector<Eigen::Vector3d> point3d_vec_;
  std::vector<Eigen::Vector3d> color_vec_;
};
}  // namespace FCA
#endif
