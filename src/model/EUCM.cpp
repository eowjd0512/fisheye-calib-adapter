#include "model/EUCM.hpp"

namespace FCA
{
namespace model
{

EUCM::EUCM(const std::string & model_name, const std::string & config_path)
: Base(model_name, config_path)
{
}

void EUCM::parse() {}
void EUCM::initialize() {}

Eigen::Vector2d EUCM::project(const Eigen::Vector3d & point3d) const
{
  const double X = point3d.x();
  const double Y = point3d.y();
  const double Z = point3d.z();

  const double d = std::sqrt(params.beta * (X * X + Y * Y) + Z * Z);
  const double denom = params.alpha * d + (1.0 - params.alpha) * Z;

  constexpr double PRECISION = 1e-3;
  if (denom < PRECISION) {
    Eigen::Vector2d point2d(-1, -1);
    return point2d;
  }

  // if (params.alpha > 0.5) {
  //   // Check that the point is in the upper hemisphere in case of ellipsoid
  //   const double zn = Z / denom;
  //   const T C = (alpha - T(1.)) / (alpha + alpha - T(1.));
  //   if (zn < C) return false;
  // }
  
  Eigen::Vector2d point2d;
  point2d.x() = params.fx * (X / denom) + params.cx;
  point2d.y() = params.fy * (Y / denom) + params.cy;
  
  return point2d;
}

Eigen::Vector3d EUCM::unproject(const Eigen::Vector2d & point2d) const
{
  const double fx = params.fx;
  const double fy = params.fy;
  const double cx = params.cx;
  const double cy = params.cy;
  const double alpha = params.alpha;
  const double beta = params.beta;
  const double u = point2d.x();
  const double v = point2d.y();

  const double mx = (u - cx) / fx;
  const double my = (v - cy) / fy;

  const double r_squared = (mx * mx) + (my * my);
  const double gamma = 1.0 - alpha;
  const double num = 1.0 - r_squared * alpha * alpha * beta;
  const double det = 1.0 - (alpha - gamma) * beta * r_squared;
  double denom = gamma + alpha * sqrt(det);
  if (det < 0) {
    denom = gamma;
  }
  const double mz = num / denom;
  const double norm = std::sqrt((mx * mx) + (my * my) + (mz * mz));

  Eigen::Vector3d point3d;
  point3d.x() = mx / norm;
  point3d.y() = my / norm;
  point3d.z() = mz / norm;

  return point3d;
}

Eigen::MatrixXd EUCM::calculate_jacobian() const {}

void EUCM::optimize() {}

void EUCM::print() const {}

}  // namespace model
}  // namespace FCA
