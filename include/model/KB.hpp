/*
MIT License

Copyright (c) 2024 Sangjun Lee, STRADVISION

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef MODEL__KB8_HPP_
#define MODEL__KB8_HPP_

#include <string>

#include "model/base.hpp"

namespace FCA
{
namespace model
{
class KB : public Base
{
public:
  struct Params
  {
    double k1;
    double k2;
    double k3;
    double k4;
  };

  KB(const std::string & model_name, const std::string & config_path);

  void parse() override;
  void set_sample_points(const std::vector<Eigen::Vector2d> & point2d_vec) override;
  void initialize(
    const Base::Params & common_params, const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec) override;
  Eigen::Vector2d project(const Eigen::Vector3d & point3d, bool condition) const override;
  Eigen::Vector3d unproject(const Eigen::Vector2d & point2d, bool condition) const override;
  void optimize(
    const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec, bool display_optimization_progress) override;
  void print() const override;
  void save_result(const std::string & result_path) const override;

  static bool check_proj_condition(double z);
  void evaluate(const model::Base * const gt) override;
  const Params & get_distortion_params() const { return distortion_; };
private:
  Params distortion_;
};

class KBAnalyticCostFunction : public ceres::SizedCostFunction<2, 8>
{
public:
  KBAnalyticCostFunction(double gt_u, double gt_v, double obs_x, double obs_y, double obs_z)
  : gt_u_(gt_u), gt_v_(gt_v), obs_x_(obs_x), obs_y_(obs_y), obs_z_(obs_z)
  {
  }

  virtual ~KBAnalyticCostFunction() {}

  virtual bool Evaluate(
    double const * const * parameters, double * residuals, double ** jacobians) const
  {
    const double fx = parameters[0][0];
    const double fy = parameters[0][1];
    const double cx = parameters[0][2];
    const double cy = parameters[0][3];
    const double k1 = parameters[0][4];
    const double k2 = parameters[0][5];
    const double k3 = parameters[0][6];
    const double k4 = parameters[0][7];

    const double r = std::sqrt(obs_x_ * obs_x_ + obs_y_ * obs_y_);
    const double theta = std::atan2(r, obs_z_);
    const double theta2 = theta * theta;
    const double theta3 = theta * theta2;
    const double theta5 = theta3 * theta2;
    const double theta7 = theta5 * theta2;
    const double theta9 = theta7 * theta2;

    const double d_theta = theta + (k1 * theta3) + (k2 * theta5) + (k3 * theta7) + (k4 * theta9);
    const double x_r = obs_x_ / r;
    const double y_r = obs_y_ / r;

    residuals[0] = (fx * d_theta * x_r + cx) - gt_u_;
    residuals[1] = (fy * d_theta * y_r + cy) - gt_v_;

    if (jacobians) {
      if (jacobians[0]) {
        Eigen::Map<Eigen::Matrix<double, 2, 8, Eigen::RowMajor>> jacobian_kb8(jacobians[0]);
        jacobian_kb8.col(0) << d_theta * x_r, 0.0;  // ∂residual_x / ∂fx, ∂residual_y / ∂fx
        jacobian_kb8.col(1) << 0.0, d_theta * y_r;  // ∂residual_x / ∂fy, ∂residual_y / ∂fy
        jacobian_kb8.col(2) << 1.0, 0.0;            // ∂residual_x / ∂cx, ∂residual_y / ∂cx
        jacobian_kb8.col(3) << 0.0, 1.0;            // ∂residual_x / ∂cy, ∂residual_y / ∂cy

        // ∂residual / ∂k1..4, ∂residual / ∂p * ∂p / ∂d(theta) * ∂d(theta) / ∂k1..4
        Eigen::Matrix2d de_dp;
        de_dp << fx, 0.0, 0.0, fy;
        Eigen::Vector2d dp_dd_theta;
        dp_dd_theta << x_r, y_r;
        Eigen::Matrix<double, 1, 4> dd_theta_dks;
        dd_theta_dks << theta3, theta5, theta7, theta9;
        jacobian_kb8.rightCols<4>() = de_dp * dp_dd_theta * dd_theta_dks;
      }
    }

    return true;
  }

private:
  const double gt_u_, gt_v_, obs_x_, obs_y_, obs_z_;
};

struct KBAutoDiffCostFunctor
{
  KBAutoDiffCostFunctor(double gt_u, double gt_v, double obs_x, double obs_y, double obs_z)
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
    T k1 = parameters[4];
    T k2 = parameters[5];
    T k3 = parameters[6];
    T k4 = parameters[7];

    T r = ceres::sqrt(T(obs_x_) * T(obs_x_) + T(obs_y_) * T(obs_y_));
    T theta = ceres::atan2(r, T(obs_z_));
    T theta2 = theta * theta;
    T theta3 = theta * theta2;
    T theta5 = theta3 * theta2;
    T theta7 = theta5 * theta2;
    T theta9 = theta7 * theta2;
    T x_r = obs_x_ / r;
    T y_r = obs_y_ / r;

    T d_theta =
      T(theta) + (k1 * T(theta3)) + (k2 * T(theta5)) + (k3 * T(theta7)) + (k4 * T(theta9));

    residuals[0] = (fx * d_theta * T(x_r) + cx) - T(gt_u_);
    residuals[1] = (fy * d_theta * T(y_r) + cy) - T(gt_v_);

    return true;
  }

private:
  const double gt_u_, gt_v_, obs_x_, obs_y_, obs_z_;
};

}  // namespace model
}  // namespace FCA
#endif
