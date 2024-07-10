/*
MIT License

Copyright (c) 2024 Sangjun Lee, STRADVISION

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

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
  struct Config
  {
    int32_t sample_point;
    bool display_optimization_progress;
  };

  Adapter(
    const Config & config, FisheyeCameraModel * const input_model,
    FisheyeCameraModel * const output_model);

  void adapt();

  void evaluate(FisheyeCameraModel * const gt_output_model = nullptr);

  cv::Mat get_image(const std::vector<Eigen::Vector2d> & point2d_vec) const;

  void display_point3d_vec(
    const std::string & name, const std::vector<Eigen::Vector3d> & point3d_vec) const;

  void set_image(const std::string & image_path);

  void show_image(const std::string & result_path);

private:
  std::vector<Eigen::Vector2d> sample_points(int32_t width, int32_t height, int n);

  void prepare_data(
    std::vector<Eigen::Vector2d> & original_point2d_vec,
    std::vector<Eigen::Vector2d> & input_point2d_vec,
    std::vector<Eigen::Vector2d> & output_point2d_vec);

  void evaluate_reprojection_error();

  double calculate_psnr(const cv::Mat & img1, const cv::Mat & img2);

  double calculate_ssim(const cv::Mat & img1, const cv::Mat & img2);

  Config config_;

  FisheyeCameraModel * const input_model_;
  FisheyeCameraModel * const output_model_;
  std::vector<Eigen::Vector2d> point2d_vec_;
  std::vector<Eigen::Vector3d> point3d_vec_;
  std::vector<cv::Vec3b> color_vec_;

  cv::Mat image_;
};
}  // namespace FCA
#endif
