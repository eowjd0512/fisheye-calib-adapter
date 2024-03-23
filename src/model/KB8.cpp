#include "model/KB8.hpp"

namespace FCA
{
namespace model
{

KB8::KB8(const std::string & model_name, const std::string & config_path)
: Base(model_name, config_path)
{
}

void KB8::parse() {}
void KB8::initialize() {}

Eigen::Vector2d KB8::project(const Eigen::Vector3d & point3d) const
{
  const double X = point3d.x();
  const double Y = point3d.y();
  const double Z = point3d.z();
  const double r = std::sqrt(X * X + Y * Y);
  const double theta = std::atan2(r, Z);
  const double theta2 = theta * theta;
  const double theta3 = theta * theta2;
  const double theta5 = theta3 * theta2;
  const double theta7 = theta5 * theta2;
  const double theta9 = theta7 * theta2;
  const double d_theta = theta + (params.k1 * theta3) + (params.k2 * theta5) +
                         (params.k3 * theta7) + (params.k4 * theta9);

  Eigen::Vector2d point2d;
  point2d.x() = params.fx * d_theta * (X / r) + params.cx;
  point2d.y() = params.fy * d_theta * (Y / r) + params.cy;

  return point2d;
}

Eigen::Vector3d KB8::unproject(const Eigen::Vector2d & point2d) const
{
  const double u = point2d.x();
  const double v = point2d.y();

  const double mx = (u - params.cx) / params.fx;
  const double my = (v - params.cy) / params.fy;

  double ru = std::sqrt(mx * mx + my * my);
  ru = std::min(std::max(-M_PI / 2.0, ru), M_PI / 2.0);
  double theta = ru;
  constexpr double PRECISION = 1e-8;
  if (ru > PRECISION) {
    // newton raphson
    while (true) {
      const double theta2 = theta * theta;
      const double theta4 = theta2 * theta2;
      const double theta6 = theta4 * theta2;
      const double theta8 = theta4 * theta4;
      const double k1_theta2 = params.k1 * theta2;
      const double k2_theta4 = params.k2 * theta4;
      const double k3_theta6 = params.k3 * theta6;
      const double k4_theta8 = params.k4 * theta8;
      const double f = theta * (1.0 + k1_theta2 + k2_theta4 + k3_theta6 + k4_theta8) - ru;
      const double f_prime =
        1.0 + (3.0 * k1_theta2) + (5.0 * k2_theta4) + (7.0 * k3_theta6) + (9.0 * k4_theta8);
      const double delta = f / f_prime;

      theta -= delta;
      if (std::abs(delta) < PRECISION) {
        break;
      }
    }

    Eigen::Vector3d point3d;
    point3d.x() = sin(theta) * (mx / ru);
    point3d.y() = sin(theta) * (my / ru);
    point3d.z() = cos(theta);

    return point3d;
  }

  Eigen::MatrixXd KB8::calculate_jacobian() const
  {
    //du/drho * drho/d[k1,k2,k3,k4] = [fx 0;0 fy] * [theta^3 * x/r] /  [theta^5 * x/r] / [theta^7 * x/r] / [theta^9 * x/r]
  }
  void KB8::optimize() {}
  void KB8::print() const {}
}  // namespace model
}  // namespace FCA
