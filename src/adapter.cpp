#include "adapter.hpp"

namespace FCA
{

Adapter::Adapter(
  const FCA::FisheyeCameraModel * const input_model, FCA::FisheyeCameraModel * const output_model)
: input_model_(input_model), output_model_(output_model)
{
}

void Adapter::adapt()
{
  const auto & common_params = input_model_->get_common_params();
  // Sample points
  constexpr auto NUM_SAMPLE_POINTS = 10000;
  const std::vector<Eigen::Vector2d> sampled_point2d_vec =
    this->sample_points(common_params.width, common_params.height, NUM_SAMPLE_POINTS);

  // Unproject samples from input model and associate <3D, 2D>
  std::vector<Eigen::Vector2d> point2d_vec;
  std::vector<Eigen::Vector3d> point3d_vec;
  point2d_vec.reserve(NUM_SAMPLE_POINTS);
  point3d_vec.reserve(NUM_SAMPLE_POINTS);
  for (const auto & point2d : sampled_point2d_vec) {
    const Eigen::Vector3d point3d = input_model_->unproject(point2d);
    if (point3d.z() > 0) {
      point2d_vec.emplace_back(point2d);
      point3d_vec.emplace_back(point3d);
    }
  }

  // Initialize output model
  output_model_->initialize(common_params, point3d_vec, point2d_vec);

  // Optimize output model
  output_model_->optimize(point3d_vec, point2d_vec);

  //TODO: draw sample points in input model

  // TODO: draw projected points output model

  double error = 0.0;
  for (auto i = 0U; i < point3d_vec.size(); ++i) {
    const auto & point3d = point3d_vec.at(i);
    const auto & point2d = point2d_vec.at(i);
    const Eigen::Vector2d projected_point2d = input_model_->project(point3d);
    // std::cout << point2d.transpose() << " " << projected_point2d.transpose() << std::endl;
    error += (point2d - projected_point2d).norm();
  }
  std::cout << error / point3d_vec.size() << std::endl;
}

void Adapter::compare_image(const std::string & image_path, bool show)
{
  // if(!output_model_->is_estimated()){
  //   this->adapt();
  // }
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
