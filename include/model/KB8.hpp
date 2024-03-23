#ifndef MODEL__EUCM_HPP_
#define MODEL__EUCM_HPP_

#include <string>

#include "model/base.hpp"

namespace FCA
{
namespace model
{
class KB8 : public Base
{
public:
  struct Params
  {
    double k1;
    double k2;
    double k3;
    double k4;
  };

  KB8(const std::string & model_name, const std::string & config_path);

  void parse() override;
  void initialize(
    const Base::Params & intrinsic, const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec) override;
  Eigen::Vector2d project(const Eigen::Vector3d & point3d) const override;
  Eigen::Vector3d unproject(const Eigen::Vector2d & point2d) const override;
  Eigen::MatrixXd calculate_jacobian() const override;
  void optimize() override;
  void print() const override;

private:
  Params distortion_;
};
}  // namespace model
}  // namespace FCA
#endif
