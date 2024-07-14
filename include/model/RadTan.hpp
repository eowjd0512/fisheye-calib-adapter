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

#ifndef MODEL__RADTAN_HPP_
#define MODEL__RADTAN_HPP_

#include <string>

#include "model/base.hpp"

namespace FCA
{
namespace model
{
class RadTan : public Base
{
public:
  struct Params
  {
    double k1;
    double k2;
    double k3;
    double p1;
    double p2;
  };

  RadTan(const std::string & model_name, const std::string & config_path);

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

class RadTanAnalyticCostFunction : public ceres::SizedCostFunction<2, 9>
{
public:
  RadTanAnalyticCostFunction(double gt_u, double gt_v, double obs_x, double obs_y, double obs_z)
  : gt_u_(gt_u), gt_v_(gt_v), obs_x_(obs_x), obs_y_(obs_y), obs_z_(obs_z)
  {
  }

  virtual ~RadTanAnalyticCostFunction() {}

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
    const double p1 = parameters[0][7];
    const double p2 = parameters[0][8];

    const double x_prime = obs_x_ / obs_z_;
    const double y_prime = obs_y_ / obs_z_;

    const double r2 = x_prime * x_prime + y_prime * y_prime;
    const double r4 = r2 * r2;
    const double r6 = r4 * r2;

    const double d = 1.0 + k1 * r2 + k2 * r4 + k3 * r6;
    const double x_distorted =
      d * x_prime + 2.0 * p1 * x_prime * y_prime + p2 * (r2 + 2.0 * x_prime * x_prime);
    const double y_distorted =
      d * y_prime + 2.0 * p2 * x_prime * y_prime + p1 * (r2 + 2.0 * y_prime * y_prime);

    residuals[0] = (fx * x_distorted + cx) - gt_u_;
    residuals[1] = (fy * y_distorted + cy) - gt_v_;

    if (jacobians) {
      if (jacobians[0]) {
        Eigen::Map<Eigen::Matrix<double, 2, 9, Eigen::RowMajor>> jacobian_radtan(jacobians[0]);
        jacobian_radtan.col(0) << x_distorted, 0.0;
        jacobian_radtan.col(1) << 0.0, y_distorted;
        jacobian_radtan.col(2) << 1.0, 0.0;
        jacobian_radtan.col(3) << 0.0, 1.0;
        jacobian_radtan.col(4) << x_prime * r2, y_prime * r2;
        jacobian_radtan.col(5) << x_prime * r4, y_prime * r4;
        jacobian_radtan.col(6) << x_prime * r6, y_prime * r6;
        jacobian_radtan.col(7) << 2.0 * x_prime * y_prime, r2 + 2.0 * y_prime * y_prime;
        jacobian_radtan.col(8) << r2 + 2.0 * x_prime * x_prime, 2.0 * x_prime * y_prime;
      }
    }

    return true;
  }

private:
  const double gt_u_, gt_v_, obs_x_, obs_y_, obs_z_;
};

}  // namespace model
}  // namespace FCA
#endif
