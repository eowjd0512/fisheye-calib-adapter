#include "adapter.hpp"

#include <nanoflann.hpp>

#include "model/base.hpp"
#include "opencv2/imgcodecs.hpp"
#include "utils.hpp"
namespace FCA
{

Adapter::Adapter(
  FCA::FisheyeCameraModel * const input_model, FCA::FisheyeCameraModel * const output_model)
: input_model_(input_model), output_model_(output_model)
{
  this->input_model_->parse();
}

void Adapter::adapt()
{
  input_model_->print();
  const auto & common_params = input_model_->get_common_params();
  // Sample points
  constexpr auto NUM_SAMPLE_POINTS = 100000;
  const std::vector<Eigen::Vector2d> sampled_point2d_vec =
    this->sample_points(common_params.width, common_params.height, NUM_SAMPLE_POINTS);

  // Unproject samples from input model and associate <3D, 2D>
  point2d_vec_.reserve(NUM_SAMPLE_POINTS);
  point3d_vec_.reserve(NUM_SAMPLE_POINTS);
  for (const auto & point2d : sampled_point2d_vec) {
    const Eigen::Vector3d point3d = input_model_->unproject(point2d);
    const auto point2d_est = input_model_->project(point3d);

    // std::cout << point3d.transpose() << std::endl;
    // std::cout << point2d.transpose() << std::endl;
    // std::cout << point2d_est.transpose() << std::endl;

    if (point3d.z() > 0.0) {
      point2d_vec_.emplace_back(point2d);
      point3d_vec_.emplace_back(point3d);
      if (!image_.empty()) {
        color_vec_.emplace_back(image_.at<cv::Vec3b>(point2d[1], point2d[0]));
      }
    }
  }

  this->display_point2d_vec("Input model's projection", point2d_vec_);
  // this->display_point3d_vec("Input model's unprojection", point3d_vec_);

  std::cout << "sampled point size: " << point3d_vec_.size() << std::endl;
  // Initialize output model
  output_model_->initialize(common_params, point3d_vec_, point2d_vec_);

  // Optimize output model
  output_model_->optimize(point3d_vec_, point2d_vec_);
  output_model_->print();

  std::vector<Eigen::Vector2d> output_point2d_vec;
  output_point2d_vec.reserve(point3d_vec_.size());
  for (const auto & point3d : point3d_vec_) {
    const Eigen::Vector2d point2d = input_model_->project(point3d);
    if (point2d.y() >= 0.0) {
      output_point2d_vec.emplace_back(point2d);
    }
  }
  this->display_point2d_vec("Output model's projection", output_point2d_vec);

  cv::Mat output_recovered_image =
    this->recover_image(common_params.width, common_params.height, point3d_vec_);
  // cv::resize(
  //   output_recovered_image, output_recovered_image,
  //   cv::Size(common_params.width / 2, common_params.height / 2));
  cv::imshow("Output recovered image", output_recovered_image);
  cv::waitKey(0);
}

void Adapter::evaluate()
{
  double error = 0.0;
  for (auto i = 0U; i < point3d_vec_.size(); ++i) {
    const auto & point3d = point3d_vec_.at(i);
    const auto & point2d = point2d_vec_.at(i);
    const Eigen::Vector2d projected_point2d = input_model_->project(point3d);
    error += (point2d - projected_point2d).norm();
  }
  std::cout << "error: " << error / point3d_vec_.size() << std::endl;
}

void Adapter::set_image(const std::string & image_path)
{
  image_ = cv::imread(image_path, cv::IMREAD_COLOR);
  if (image_.empty()) {
    ERROR_STR("Image read failed");
  }
}

void Adapter::display_point2d_vec(
  const std::string & name, const std::vector<Eigen::Vector2d> & point2d_vec) const
{
  const double width = input_model_->get_common_params().width;
  const double height = input_model_->get_common_params().height;

  cv::Mat image(height, width, CV_8UC3, cv::Scalar(0, 0, 0));

  for (const auto & point : point2d_vec) {
    cv::Point cvPoint(point[0], point[1]);
    cv::Scalar color = cv::Scalar(255, 255, 255);
    if (!image_.empty()) {
      color = image_.at<cv::Vec3b>(point[1], point[0]);
    }
    cv::circle(image, cvPoint, 3, color, -1);
  }

  // cv::resize(image, image, cv::Size(width / 2, height / 2));
  cv::imshow(name, image);
  cv::waitKey(0);
}

void Adapter::display_point3d_vec(
  const std::string & name, const std::vector<Eigen::Vector3d> & point3d_vec) const
{
  std::vector<cv::Point3f> cv_points_3d;
  std::vector<cv::Vec3b> cv_colors;
  cv_points_3d.reserve(point3d_vec_.size());
  cv_colors.reserve(point3d_vec_.size());
  for (auto i = 0U; i < point3d_vec_.size(); ++i) {
    const auto & point3d = point3d_vec_.at(i);

    cv_points_3d.emplace_back(point3d.x(), point3d.y(), point3d.z());

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
