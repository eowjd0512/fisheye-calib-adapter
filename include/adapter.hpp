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
    const Config & config, FCA::FisheyeCameraModel * const input_model,
    FCA::FisheyeCameraModel * const output_model);

  void adapt();

  void evaluate();

  cv::Mat get_image(const std::vector<Eigen::Vector2d> & point2d_vec) const;

  void display_point3d_vec(
    const std::string & name, const std::vector<Eigen::Vector3d> & point3d_vec) const;

  cv::Mat recover_image(
    double width, double height, const std::vector<Eigen::Vector3d> & point3d_vec) const;

  void set_image(const std::string & image_path);

  void show_image();

private:
  std::vector<Eigen::Vector2d> sample_points(int32_t width, int32_t height, int n);

  void prepare_data(
    std::vector<Eigen::Vector2d> & original_point2d_vec,
    std::vector<Eigen::Vector2d> & input_point2d_vec,
    std::vector<Eigen::Vector2d> & output_point2d_vec);

  void evaluate_reprojection_error();

  void evaluate_psnr(
    const cv::Mat & original_img, const cv::Mat & input_model_img,
    const cv::Mat & output_model_img);

  void evaluate_ssim(
    const cv::Mat & original_img, const cv::Mat & input_model_img,
    const cv::Mat & output_model_img);

  double calculate_psnr(const cv::Mat & img1, const cv::Mat & img2);

  double calculate_ssim(const cv::Mat & img1, const cv::Mat & img2);

  Config config_;

  FCA::FisheyeCameraModel * const input_model_;
  FCA::FisheyeCameraModel * const output_model_;
  std::vector<Eigen::Vector2d> point2d_vec_;
  std::vector<Eigen::Vector3d> point3d_vec_;
  std::vector<cv::Vec3b> color_vec_;

  cv::Mat image_;
};
}  // namespace FCA
#endif
