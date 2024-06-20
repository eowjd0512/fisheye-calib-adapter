#ifndef MODEL__EUCM_HPP_
#define MODEL__EUCM_HPP_

#include <string>

#include "model/base.hpp"

namespace FCA
{
namespace model
{
class EUCM : public Base
{
public:
  struct Params
  {
    double alpha;
    double beta;
  };

  EUCM(const std::string & model_name, const std::string & config_path);

  void parse() override;
  void set_sample_points(const std::vector<Eigen::Vector2d> & point2d_vec) override;
  void initialize(
    const Base::Params & common_params, const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec) override;
  Eigen::Vector2d project(const Eigen::Vector3d & point3d) const override;
  Eigen::Vector3d unproject(const Eigen::Vector2d & point2d) const override;
  void optimize(
    const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec) override;
  void print() const override;
  void save_result(const std::string & result_path) const override;

  static bool check_proj_condition(double z, double denom, double alpha, double beta);
  static bool check_unproj_condition(double r_squared, double alpha, double beta);

private:
  Params distortion_;
};

class EUCMAnalyticCostFunction : public ceres::SizedCostFunction<2, 6>
{
public:
  EUCMAnalyticCostFunction(double gt_u, double gt_v, double obs_x, double obs_y, double obs_z)
  : gt_u_(gt_u), gt_v_(gt_v), obs_x_(obs_x), obs_y_(obs_y), obs_z_(obs_z)
  {
  }

  virtual ~EUCMAnalyticCostFunction() {}

  virtual bool Evaluate(
    double const * const * parameters, double * residuals, double ** jacobians) const
  {
    const double fx = parameters[0][0];
    const double fy = parameters[0][1];
    const double cx = parameters[0][2];
    const double cy = parameters[0][3];
    const double alpha = parameters[0][4];
    const double beta = parameters[0][5];

    const double u_cx = gt_u_ - cx;
    const double v_cy = gt_v_ - cy;
    const double r_squared = obs_x_ * obs_x_ + obs_y_ * obs_y_;
    const double d = std::sqrt(beta * (r_squared) + obs_z_ * obs_z_);
    const double denom = alpha * d + (1.0 - alpha) * obs_z_;

    constexpr double PRECISION = 1e-3;
    if ((denom < PRECISION) || !EUCM::check_proj_condition(obs_z_, denom, alpha, beta)) {
      return false;
    }

    residuals[0] = fx * obs_x_ - u_cx * denom;
    residuals[1] = fy * obs_y_ - v_cy * denom;

    if (jacobians) {
      if (jacobians[0]) {
        Eigen::Map<Eigen::Matrix<double, 2, 6, Eigen::RowMajor>> jacobian_eucm(jacobians[0]);
        jacobian_eucm.col(0) << obs_x_, 0.0;  // ∂residual_x / ∂fx, ∂residual_y / ∂fx
        jacobian_eucm.col(1) << 0.0, obs_y_;  // ∂residual_x / ∂fy, ∂residual_y / ∂fy
        jacobian_eucm.col(2) << denom, 0.0;   // ∂residual_x / ∂cx, ∂residual_y / ∂cx
        jacobian_eucm.col(3) << 0.0, denom;   // ∂residual_x / ∂cy, ∂residual_y / ∂cy
        jacobian_eucm.col(4) << (obs_z_ - d) * u_cx,
          (obs_z_ - d) * v_cy;  // ∂residual_x / ∂alpha, ∂residual_y / ∂alpha
        jacobian_eucm.col(5) << -(alpha * r_squared * u_cx) / (2.0 * d),
          -(alpha * r_squared * v_cy) / (2.0 * d);  // ∂residual_x / ∂beta, ∂residual_y / ∂beta
      }
    }

    return true;
  }

private:
  const double gt_u_, gt_v_, obs_x_, obs_y_, obs_z_;
};

struct EUCMAutoDiffCostFunctor
{
  EUCMAutoDiffCostFunctor(double gt_u, double gt_v, double obs_x, double obs_y, double obs_z)
  : gt_u_(gt_u), gt_v_(gt_v), obs_x_(obs_x), obs_y_(obs_y), obs_z_(obs_z)
  {
  }

  template <typename T>
  bool operator()(const T * const parameters, T * residuals) const
  {
    T fx = parameters[0];
    T fy = parameters[1];
    T cx = parameters[2];
    T cy = parameters[3];
    T alpha = parameters[4];
    T beta = parameters[5];

    const double r_squared = obs_x_ * obs_x_ + obs_y_ * obs_y_;
    T u_cx = gt_u_ - cx;
    T v_cy = gt_v_ - cy;
    T d = ceres::sqrt(beta * T(r_squared) + T(obs_z_) * T(obs_z_));
    T denom = alpha * d + (1.0 - alpha) * T(obs_z_);

    constexpr double PRECISION = 1e-3;
    if ((denom < PRECISION) || !EUCM::check_proj_condition(obs_z_, denom, alpha, beta)) {
      return false;
    }

    residuals[0] = fx * T(obs_x_) - u_cx * denom;
    residuals[1] = fy * T(obs_y_) - v_cy * denom;

    return true;
  }

private:
  const double gt_u_, gt_v_, obs_x_, obs_y_, obs_z_;
};

}  // namespace model
}  // namespace FCA
#endif
