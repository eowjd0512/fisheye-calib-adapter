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

#include "adapter.hpp"

#include "model/base.hpp"
#include "opencv2/imgcodecs.hpp"
#include "utils.hpp"

namespace FCA
{

Adapter::Adapter(
  const Config & config, FisheyeCameraModel * const input_model,
  FisheyeCameraModel * const output_model)
: config_(config), input_model_(input_model), output_model_(output_model)
{
  this->input_model_->parse();
}

void Adapter::adapt()
{
  input_model_->print();
  const auto & common_params = input_model_->get_common_params();

  // Sample points
  const std::vector<Eigen::Vector2d> sampled_point2d_vec =
    this->sample_points(common_params.width, common_params.height, config_.sample_point);

  // Unproject samples from input model and associate <3D, 2D>
  point2d_vec_.reserve(config_.sample_point);
  point3d_vec_.reserve(config_.sample_point);
  for (const auto & point2d : sampled_point2d_vec) {
    const Eigen::Vector3d point3d = input_model_->unproject(point2d);

    if (point3d.z() > 0.0) {
      point2d_vec_.emplace_back(point2d);
      point3d_vec_.emplace_back(point3d);
      if (!image_.empty()) {
        color_vec_.emplace_back(image_.at<cv::Vec3b>(point2d[1], point2d[0]));
      }
    }
  }

  // Initialize output model
  output_model_->set_sample_points(sampled_point2d_vec);  // only for OcamCalib
  output_model_->initialize(common_params, point3d_vec_, point2d_vec_);

  // Optimize output model
  output_model_->optimize(point3d_vec_, point2d_vec_, config_.display_optimization_progress);
  output_model_->print();

  // display_point3d_vec("3d points", point3d_vec_);
}

void Adapter::show_image(const std::string & result_path)
{
  std::vector<Eigen::Vector2d> original_point2d_vec;
  std::vector<Eigen::Vector2d> input_point2d_vec;
  std::vector<Eigen::Vector2d> input_to_output_point2d_vec;

  this->prepare_data(original_point2d_vec, input_point2d_vec, input_to_output_point2d_vec);

  cv::Mat original_image = this->get_image(original_point2d_vec);
  cv::Mat input_model_image = this->get_image(input_point2d_vec);
  cv::Mat input_output_model_image = this->get_image(input_to_output_point2d_vec);

  cv::imshow("Original image", original_image);
  cv::imshow("Input model's projection", input_model_image);
  cv::imshow("Output model's projection", input_output_model_image);

  cv::imwrite(result_path + "/original_image.png", original_image);
  cv::imwrite(result_path + "/input_model_image.png", input_model_image);
  cv::imwrite(result_path + "/input_to_output_model_image.png", input_output_model_image);
  cv::waitKey(0);
}
void Adapter::prepare_data(
  std::vector<Eigen::Vector2d> & original_point2d_vec,
  std::vector<Eigen::Vector2d> & input_point2d_vec,
  std::vector<Eigen::Vector2d> & output_point2d_vec)
{
  const auto & common_params = output_model_->get_common_params();
  for (auto x = 0; x < common_params.width; ++x) {
    for (auto y = 0; y < common_params.height; ++y) {
      Eigen::Vector2d point2d(x, y);
      const Eigen::Vector3d input_point3d = input_model_->unproject(point2d, false);
      // const Eigen::Vector3d output_point3d = output_model_->unproject(point2d, false);
      const Eigen::Vector2d input_point2d = input_model_->project(input_point3d, false);
      const Eigen::Vector2d output_point2d = output_model_->project(input_point3d, false);
      original_point2d_vec.emplace_back(point2d);
      input_point2d_vec.emplace_back(input_point2d);
      output_point2d_vec.emplace_back(output_point2d);
    }
  }
}

void Adapter::evaluate(FisheyeCameraModel * const gt_output_model)
{
  evaluate_reprojection_error();

  if (!image_.empty()) {
    std::vector<Eigen::Vector2d> original_point2d_vec;
    std::vector<Eigen::Vector2d> input_point2d_vec;
    std::vector<Eigen::Vector2d> output_point2d_vec;
    this->prepare_data(original_point2d_vec, input_point2d_vec, output_point2d_vec);

    cv::Mat original_img = this->get_image(original_point2d_vec);
    cv::Mat input_model_img = this->get_image(input_point2d_vec);
    cv::Mat input_to_output_model_img = this->get_image(output_point2d_vec);

    std::cout << "psnr from input model to output model: "
              << calculate_psnr(original_img, input_to_output_model_img) << std::endl;
    std::cout << "ssim from input model to output model: "
              << calculate_ssim(original_img, input_to_output_model_img) << std::endl;
  }

  if (gt_output_model) {
    gt_output_model->parse();
    output_model_->evaluate(gt_output_model);
  }
}

void Adapter::evaluate_reprojection_error()
{
  double output_model_error = 0.0;

  for (auto i = 0U; i < point3d_vec_.size(); ++i) {
    const auto & point3d = point3d_vec_.at(i);
    const auto & point2d = point2d_vec_.at(i);
    const Eigen::Vector2d output_point2d = output_model_->project(point3d);
    output_model_error += (point2d - output_point2d).norm();
  }
  std::cout << "reprojection error from input model to output model: "
            << output_model_error / point3d_vec_.size() << std::endl;
}

double Adapter::calculate_psnr(const cv::Mat & img1, const cv::Mat & img2)
{
  cv::Mat gray1, gray2;
  cv::cvtColor(img1, gray1, cv::COLOR_BGR2GRAY);
  cv::cvtColor(img2, gray2, cv::COLOR_BGR2GRAY);

  cv::Mat mask = (gray1 != 0) & (gray2 != 0);

  cv::Mat s1;

  cv::absdiff(img1, img2, s1);
  s1.convertTo(s1, CV_32F);

  cv::Mat s1_masked;
  s1.copyTo(s1_masked, mask);

  cv::Mat s1_squared = s1_masked.mul(s1_masked);

  cv::Scalar s = cv::sum(s1_squared);
  double sse = s.val[0];

  if (sse <= 1e-10) {
    return 0;
  } else {
    double mse = sse / cv::countNonZero(mask);
    double psnr = 10.0 * log10((255 * 255) / mse);
    return psnr;
  }
}

double Adapter::calculate_ssim(const cv::Mat & img1, const cv::Mat & img2)
{
  const double C1 = 6.5025, C2 = 58.5225;
  cv::Mat mask = (img1 != 0) & (img2 != 0);

    cv::Mat I1, I2;
    img1.convertTo(I1, CV_32F);
    img2.convertTo(I2, CV_32F);

    cv::Mat I1_masked, I2_masked;
    I1.copyTo(I1_masked, mask);
    I2.copyTo(I2_masked, mask);

    cv::Mat I1_2 = I1_masked.mul(I1_masked);
    cv::Mat I2_2 = I2_masked.mul(I2_masked);
    cv::Mat I1_I2 = I1_masked.mul(I2_masked);

    cv::Mat mu1, mu2;
    cv::GaussianBlur(I1_masked, mu1, cv::Size(11, 11), 1.5);
    cv::GaussianBlur(I2_masked, mu2, cv::Size(11, 11), 1.5);

    cv::Mat mu1_2 = mu1.mul(mu1);
    cv::Mat mu2_2 = mu2.mul(mu2);
    cv::Mat mu1_mu2 = mu1.mul(mu2);

    cv::Mat sigma1_2, sigma2_2, sigma12;

    cv::GaussianBlur(I1_2, sigma1_2, cv::Size(11, 11), 1.5);
    sigma1_2 -= mu1_2;

    cv::GaussianBlur(I2_2, sigma2_2, cv::Size(11, 11), 1.5);
    sigma2_2 -= mu2_2;

    cv::GaussianBlur(I1_I2, sigma12, cv::Size(11, 11), 1.5);
    sigma12 -= mu1_mu2;

    cv::Mat t1, t2, t3;
    t1 = 2.0 * mu1_mu2 + C1;
    t2 = 2.0 * sigma12 + C2;
    t3 = t1.mul(t2);

    t1 = mu1_2 + mu2_2 + C1;
    t2 = sigma1_2 + sigma2_2 + C2;
    t1 = t1.mul(t2);

    cv::Mat ssim_map;
    cv::divide(t3, t1, ssim_map);

    cv::Scalar mssim = cv::mean(ssim_map);
    return mssim.val[0];
}

void Adapter::set_image(const std::string & image_path)
{
  image_ = cv::imread(image_path, cv::IMREAD_COLOR);
  if (image_.empty()) {
    std::cerr << "Image read failed" << std::endl;
  }
}

cv::Mat Adapter::get_image(const std::vector<Eigen::Vector2d> & point2d_vec) const
{
  const double width = input_model_->get_common_params().width;
  const double height = input_model_->get_common_params().height;

  cv::Mat image(height, width, CV_8UC3, cv::Scalar(0, 0, 0));

  for (const auto & point : point2d_vec) {
    const auto u = std::floor(point[0] + 0.5);
    const auto v = std::floor(point[1] + 0.5);

    if ((u >= 0.0) && (u < width) && (v >= 0.0) && (v < height)) {
      cv::Point cvPoint(point[0], point[1]);
      cv::Scalar color = cv::Scalar(255, 255, 255);
      if (!image_.empty()) {
        color = image_.at<cv::Vec3b>(point[1], point[0]);
      }
      cv::circle(image, cvPoint, 1, color, -1);
    }
  }
  return image;
}

void Adapter::display_point3d_vec(
  const std::string & name, const std::vector<Eigen::Vector3d> & point3d_vec) const
{
  std::vector<cv::Point3f> cv_points_3d;
  std::vector<cv::Vec3b> cv_colors;
  cv_points_3d.reserve(point3d_vec_.size());
  cv_colors.reserve(point3d_vec_.size());
  constexpr auto SCALE = 3.0;
  for (auto i = 0U; i < point3d_vec_.size(); ++i) {
    const auto & point3d = point3d_vec_.at(i);

    cv_points_3d.emplace_back(point3d.x() * SCALE, point3d.y() * SCALE, point3d.z() * SCALE);

    if (!color_vec_.empty()) {
      const auto & color = color_vec_.at(i);
      cv_colors.emplace_back(color);
    } else {
      cv_colors.emplace_back(255, 255, 255);
    }
  }

  cv::viz::Viz3d window(name);

  cv::viz::WCloud cloud_widget(cv_points_3d, cv_colors);

  window.showWidget("Coordinate", cv::viz::WCoordinateSystem());
  window.showWidget("3d points", cloud_widget);
  window.spin();
}

std::vector<Eigen::Vector2d> Adapter::sample_points(int32_t width, int32_t height, int n)
{
  std::vector<Eigen::Vector2d> points;
  const int32_t num_cells_x = std::round(std::sqrt(n * (width / static_cast<float>(height))));
  const int32_t num_cells_y = std::round(std::sqrt(n * (height / static_cast<float>(width))));

  const double cell_width = static_cast<double>(width) / static_cast<double>(num_cells_x);
  const double cell_height = static_cast<double>(height) / static_cast<double>(num_cells_y);

  for (auto i = 0; i < num_cells_y; ++i) {
    for (auto j = 0; j < num_cells_x; ++j) {
      Eigen::Vector2d p;
      p.x() = (j + 0.5f) * cell_width;
      p.y() = (i + 0.5f) * cell_height;
      points.push_back(p);
    }
  }

  return points;
}
}  // namespace FCA
