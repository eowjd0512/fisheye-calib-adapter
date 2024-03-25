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

void KB8::initialize(
  const Base::Params & intrinsic, const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
  assert(point3d_vec.size() == point2d_vec.size());

  // set fx,fy,cx,cy
  intrinsic_ = intrinsic;

  // set k1, k2, k3, k4
  Eigen::MatrixXd A(point3d_vec.size(), 4);
  Eigen::VectorXd b(point3d_vec.size(), 1);

  for (auto i = 0U; i < point3d_vec.size(); ++i) {
    const auto & point3d = point3d_vec.at(i);
    const auto & point2d = point2d_vec.at(i);
    const double X = point3d.x();
    const double Y = point3d.y();
    const double Z = point3d.z();
    const double u = point2d.x();
    const double v = point2d.y();
    const double r = std::sqrt((X * X) + (Y * Y));
    const double theta = std::atan2(r, Z);
    const double theta2 = theta * theta;
    const double theta3 = theta * theta2;
    const double theta5 = theta3 * theta2;
    const double theta7 = theta5 * theta2;
    const double theta9 = theta7 * theta2;

    A(i, 0) = theta3;
    A(i, 1) = theta5;
    A(i, 2) = theta7;
    A(i, 3) = theta9;
    b[i] = (u - intrinsic_.cx) * (r / (intrinsic_.fx * X)) - theta;
  }

  const Eigen::VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
  distortion_.k1 = x[0];
  distortion_.k2 = x[1];
  distortion_.k3 = x[2];
  distortion_.k4 = x[3];
}

Eigen::Vector2d KB8::project(const Eigen::Vector3d & point3d) const
{
  const double X = point3d.x();
  const double Y = point3d.y();
  const double Z = point3d.z();
  const double r = std::sqrt((X * X) + (Y * Y));
  const double theta = std::atan2(r, Z);
  const double theta2 = theta * theta;
  const double theta3 = theta * theta2;
  const double theta5 = theta3 * theta2;
  const double theta7 = theta5 * theta2;
  const double theta9 = theta7 * theta2;
  const double d_theta = theta + (distortion_.k1 * theta3) + (distortion_.k2 * theta5) +
                         (distortion_.k3 * theta7) + (distortion_.k4 * theta9);

  Eigen::Vector2d point2d;
  point2d.x() = intrinsic_.fx * d_theta * (X / r) + intrinsic_.cx;
  point2d.y() = intrinsic_.fy * d_theta * (Y / r) + intrinsic_.cy;

  return point2d;
}

Eigen::Vector3d KB8::unproject(const Eigen::Vector2d & point2d) const
{
  const double u = point2d.x();
  const double v = point2d.y();

  const double mx = (u - intrinsic_.cx) / intrinsic_.fx;
  const double my = (v - intrinsic_.cy) / intrinsic_.fy;

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
      const double k1_theta2 = distortion_.k1 * theta2;
      const double k2_theta4 = distortion_.k2 * theta4;
      const double k3_theta6 = distortion_.k3 * theta6;
      const double k4_theta8 = distortion_.k4 * theta8;
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

  void KB8::optimize() {}
  void KB8::print() const {}
}  // namespace model
}  // namespace FCA
