#include "adapter.hpp"

#include <nanoflann.hpp>

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
  output_model_->set_sample_points(sampled_point2d_vec);  // only for OcamLib
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
  return cv::PSNR(img1, img2);
}

double Adapter::calculate_ssim(const cv::Mat & img1, const cv::Mat & img2)
{
  const double C1 = 6.5025, C2 = 58.5225;
  int d = CV_32F;

  cv::Mat I1, I2;
  img1.convertTo(I1, d);
  img2.convertTo(I2, d);

  cv::Mat I1_2 = I1.mul(I1);   // I1^2
  cv::Mat I2_2 = I2.mul(I2);   // I2^2
  cv::Mat I1_I2 = I1.mul(I2);  // I1 * I2

  cv::Mat mu1, mu2;
  cv::GaussianBlur(I1, mu1, cv::Size(11, 11), 1.5);
  cv::GaussianBlur(I2, mu2, cv::Size(11, 11), 1.5);

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
  t3 = t1.mul(t2);  // t3 = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))

  t1 = mu1_2 + mu2_2 + C1;
  t2 = sigma1_2 + sigma2_2 + C2;
  t1 = t1.mul(t2);  // t1 = ((mu1_2 + mu2_2 + C1).*(sigma1_2 + sigma2_2 + C2))

  cv::Mat ssim_map;
  cv::divide(t3, t1, ssim_map);  // ssim_map = t3./t1;
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
    const auto u = point[0];
    const auto v = point[1];

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

cv::Mat Adapter::recover_image(
  double width, double height, const std::vector<Eigen::Vector3d> & point3d_vec) const
{
  cv::Mat output_image(height, width, CV_8UC3, cv::Scalar(0, 0, 0));

  PointCloud<double> cloud;
  cloud.pts.reserve(point3d_vec.size());
  for (const auto & pt : point3d_vec) {
    cloud.pts.push_back({pt.x(), pt.y(), pt.z()});
  }

  using MyKDTreeType = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>, PointCloud<double>, 3 /* dim */
    >;

  std::cout << "Begin to constructing tree" << std::endl;
  // register point3d_vec to flann tree
  MyKDTreeType index(
    3 /*dim*/, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
  index.buildIndex();
  std::cout << "End to constructing tree" << std::endl;

  for (auto y = 0; y < height; ++y) {
    for (auto x = 0; x < width; ++x) {
      Eigen::Vector2d point2d{x, y};
      const Eigen::Vector3d point3d = output_model_->unproject(point2d);
      if (point3d.z() > 0.0) {
        // find closest point from tree
        double query_pt[3] = {point3d.x(), point3d.y(), point3d.z()};
        size_t ret_index;
        double out_dist_sqr;
        nanoflann::KNNResultSet<double> resultSet(1);
        resultSet.init(&ret_index, &out_dist_sqr);
        index.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

        // if the distance is greater than the threshold, reject
        // else get the index of found point
        constexpr auto THRESHOLD = 1e-5;
        if (out_dist_sqr > THRESHOLD) continue;

        // get the value of color from color_vec_ using the found index
        cv::Scalar color(255., 255., 255.);
        if (!color_vec_.empty()) {
          color = color_vec_[ret_index];
        }
        // assign the color at the (y, x) position of the output_image
        output_image.at<cv::Vec3b>(y, x) = cv::Vec3b(color[0], color[1], color[2]);
      }
    }
  }
  return output_image;
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
