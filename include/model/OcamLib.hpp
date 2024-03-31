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
    const Base::Params & common_params, const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec) override;
  Eigen::Vector2d project(const Eigen::Vector3d & point3d) const override;
  Eigen::Vector3d unproject(const Eigen::Vector2d & point2d) const override;
  void optimize(
    const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec) override;
  void print() const override;
  void save_result(const std::string & result_path) const override;

  void estimate_projection_coefficients(
    const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec);

  double calculate_average_error(
    const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec);

  static bool check_proj_condition(double z);

private:
  Params distortion_;
};

class OcamLibAnalyticCostFunction : public ceres::SizedCostFunction<2, 7>
{
public:
  OcamLibAnalyticCostFunction(double gt_u, double gt_v, double obs_x, double obs_y, double obs_z)
  : gt_u_(gt_u), gt_v_(gt_v), obs_x_(obs_x), obs_y_(obs_y), obs_z_(obs_z)
  {
  }

  virtual ~OcamLibAnalyticCostFunction() {}

  virtual bool Evaluate(
    double const * const * parameters, double * residuals, double ** jacobians) const
  {
    const double c = parameters[0][0];
    const double d = parameters[0][1];
    const double e = parameters[0][2];
    const double cx = parameters[0][3];
    const double cy = parameters[0][4];
    const double k0 = parameters[0][5];
    const double k1 = parameters[0][6];
    const double k2 = parameters[0][7];
    const double k3 = parameters[0][8];
    const double k4 = parameters[0][9];

    const double u_cx = gt_u_ - cx;
    const double v_cy = gt_v_ - cy;
    const double px = obs_x_ / obs_z_;
    const double py = obs_y_ / obs_z_;
    const double r = std::sqrt(u_cx * u_cx + v_cy * v_cy);
    const double r2 = r * r;
    const double r3 = r * r2;
    const double r4 = r2 * r2;

    const double theta = k0 + k1 * r + k2 * r2 + k3 * r3 + k4 * r4;

    const double theta_px = theta * px;
    const double theta_py = theta * py;
    residuals[0] = u_cx - theta_px;
    residuals[1] = v_cy - theta_py;

    if (jacobians) {
      if (jacobians[0]) {
        Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_ocamlib(jacobians[0]);
        jacobian_ocamlib.col(0) << -1.0, 0.0;  // ∂residual_x / ∂cx, ∂residual_y / ∂cx
        jacobian_ocamlib.col(1) << 0.0, -1.0;  // ∂residual_x / ∂cy, ∂residual_y / ∂cy
        // ∂residual / ∂k0..4, ∂residual / ∂theta * ∂theta / ∂k0..4
        Eigen::Vector2d de_dtheta;
        de_dtheta << -px, -py;
        Eigen::Matrix<double, 1, 5> dtheta_dks;
        dtheta_dks << 1, r, r2, r3, r4;
        jacobian_ocamlib.rightCols<5>() = de_dtheta * dtheta_dks;
      }
    }

    return true;
  }

private:
  const double gt_u_, gt_v_, obs_x_, obs_y_, obs_z_;
};

struct OcamLibAutoDiffCostFunctor
{
  OcamLibAutoDiffCostFunctor(double gt_u, double gt_v, double obs_x, double obs_y, double obs_z)
  : gt_u_(gt_u), gt_v_(gt_v), obs_x_(obs_x), obs_y_(obs_y), obs_z_(obs_z)
  {
  }

  template <typename T>
  bool operator()(const T * const parameters, T * residuals) const
  {
    T c = parameters[0];
    T d = parameters[1];
    T e = parameters[2];
    T cx = parameters[3];
    T cy = parameters[4];
    T k0 = parameters[5];
    T k1 = parameters[6];
    T k2 = parameters[7];
    T k3 = parameters[8];
    T k4 = parameters[9];

    T c_de = c - (d * e);
    T u_cx = T(gt_u_) - cx;
    T v_cy = T(gt_v_) - cy;
    T px = T(obs_x_) / T(obs_z_);
    T py = T(obs_y_) / T(obs_z_);

    T a = u_cx - (d * v_cy);
    T b = (-e * u_cx) + (c * v_cy);
    T mx = a / c_de;
    T my = b / c_de;
    T r = ceres::sqrt(mx * mx + my * my);
    T r2 = r * r;
    T r3 = r * r2;
    T r4 = r2 * r2;

    T theta = k0 + k1 * r + k2 * r2 + k3 * r3 + k4 * r4;

    T theta_px = theta * T(px);
    T theta_py = theta * T(py);

    residuals[0] = a - theta_px * c_de;
    residuals[1] = b - theta_py * c_de;

    return true;
  }

private:
  const double gt_u_, gt_v_, obs_x_, obs_y_, obs_z_;
};

}  // namespace model
}  // namespace FCA
#endif
