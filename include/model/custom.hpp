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
  void initialize(
    const Base::Params & common_params, const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec) override;
  Eigen::Vector2d project(const Eigen::Vector3d & point3d) const override;
  Eigen::Vector3d unproject(const Eigen::Vector2d & point2d) const override;
  void optimize(
    const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec) override;
  void print() const override;

private:
  Params distortion_;
};

class CustomAnalyticCostFunction : public ceres::SizedCostFunction<2, 6>
{
public:
  CustomAnalyticCostFunction(double gt_u, double gt_v, double obs_x, double obs_y, double obs_z)
  : gt_u_(gt_u), gt_v_(gt_v), obs_x_(obs_x), obs_y_(obs_y), obs_z_(obs_z)
  {
  }

  virtual ~CustomAnalyticCostFunction() {}

  virtual bool Evaluate(
    double const * const * parameters, double * residuals, double ** jacobians) const
  {
    return true;
  }

private:
  const double gt_u_, gt_v_, obs_x_, obs_y_, obs_z_;
};

struct CustomAutoDiffCostFunctor
{
  CustomAutoDiffCostFunctor(double gt_u, double gt_v, double obs_x, double obs_y, double obs_z)
  : gt_u_(gt_u), gt_v_(gt_v), obs_x_(obs_x), obs_y_(obs_y), obs_z_(obs_z)
  {
  }

  template <typename T>
  bool operator()(const T * const parameters, T * residuals) const
  {
    return true;
  }

private:
  const double gt_u_, gt_v_, obs_x_, obs_y_, obs_z_;
};
}  // namespace model
}  // namespace FCA
#endif
