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
void Custom::initialize(
  const Base::Params & intrinsic, const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
    assert(point3d_vec.size() == point2d_vec.size());
}
Eigen::Vector2d Custom::project(const Eigen::Vector3d & point3d) const { return Eigen::Vector2d(); }
Eigen::Vector3d Custom::unproject(const Eigen::Vector2d & point2d) const
{
  return Eigen::Vector3d();
}
Eigen::MatrixXd Custom::calculate_jacobian() const {}
void Custom::optimize() {}
void Custom::print() const {}
}  // namespace model
}  // namespace FCA
