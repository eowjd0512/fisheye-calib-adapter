#ifndef MODEL__OCAMLIB_HPP_
#define MODEL__OCAMLIB_HPP_

#include <string>

#include "model/base.hpp"

namespace FCA
{
namespace model
{
class OcamLib : public Base
{
public:
  struct Params
  {
    double c;
    double d;
    double e;
    std::vector<double> proj_coeffs;
    std::vector<double> unproj_coeffs;
  };

  OcamLib(const std::string & model_name, const std::string & config_path);

  void parse() override;
  void initialize(
    const Base::Params & intrinsic, const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec) override;
  Eigen::Vector2d project(const Eigen::Vector3d & point3d) const override;
  Eigen::Vector3d unproject(const Eigen::Vector2d & point2d) const override;
  Eigen::MatrixXd calculate_jacobian() const override;
  void optimize() override;
  void print() const override;

  void estimate_projection_coefficients(
    const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec);

  // TODO: define distances for evaluation
  double evaluate(
    const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec);

private:

  Params distortion_;
};
}  // namespace model
}  // namespace FCA
#endif
