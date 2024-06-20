#include "model/custom.hpp"

namespace FCA
{
namespace model
{

Custom::Custom(const std::string & model_name, const std::string & config_path)
: Base(model_name, config_path)
{
}

void Custom::parse() {}
void Custom::set_sample_points(const std::vector<Eigen::Vector2d> & point2d_vec){};
void Custom::initialize(
  const Base::Params & common_params, const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
  assert(point3d_vec.size() == point2d_vec.size());
}
Eigen::Vector2d Custom::project(const Eigen::Vector3d & point3d, bool condition) const { return Eigen::Vector2d(); }
Eigen::Vector3d Custom::unproject(const Eigen::Vector2d & point2d, bool condition) const
{
  return Eigen::Vector3d();
}
void Custom::optimize(
  const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec,
  bool display_optimization_progress)
{
}
void Custom::print() const {}

void Custom::save_result(const std::string& result_path) const {}
}  // namespace model
}  // namespace FCA
