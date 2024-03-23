#ifndef MODEL__CUSTOM_HPP_
#define MODEL__CUSTOM_HPP_

#include <string>

#include "model/base.hpp"

namespace FCA
{
namespace model
{
class Custom : public Base
{
public:
  struct Params
  {
  };

  Custom(const std::string & model_name, const std::string & config_path);

  void parse() override;
  void initialize() override;
  Eigen::Vector2d project(const Eigen::Vector3d & point3d) const override;
  Eigen::Vector3d unproject(const Eigen::Vector2d & point2d) const override;
  Eigen::MatrixXd calculate_jacobian() const override;
  void optimize() override;
  void print() const override;

private:
  Params params;
};
}  // namespace model
}  // namespace FCA
#endif
