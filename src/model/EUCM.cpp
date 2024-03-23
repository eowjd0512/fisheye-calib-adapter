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
void EUCM::initialize(
  const Base::Params & intrinsic, const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
  assert(point3d_vec.size() == point2d_vec.size());

  // set fx,fy,cx,cy
  intrinsic_ = intrinsic;

  // set beta
  distortion_.beta = 1.0;

  // set alpha
  const int32_t mid_index = point3d_vec.size() / 2;
  const auto & point3d_mid = point3d_vec.at(mid_index);
  const auto & point2d_mid = point2d_vec.at(mid_index);

  const double X = point3d_mid.x();
  const double Y = point3d_mid.y();
  const double Z = point3d_mid.z();
  const double u = point2d_mid.x();
  const double v = point2d_mid.y();

  const double d = std::sqrt((X * X) + (Y * Y) + (Z * Z));
  const double pc = u - intrinsic_.cx;
  distortion_.alpha = ((intrinsic_.fx * X) - (pc * Z)) / (pc * (d - Z));
}

Eigen::Vector2d EUCM::project(const Eigen::Vector3d & point3d) const
{
  const double X = point3d.x();
  const double Y = point3d.y();
  const double Z = point3d.z();

  const double d = std::sqrt(distortion_.beta * ((X * X) + (Y * Y)) + (Z * Z));
  const double denom = distortion_.alpha * d + (1.0 - distortion_.alpha) * Z;

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
  point2d.x() = intrinsic_.fx * (X / denom) + intrinsic_.cx;
  point2d.y() = intrinsic_.fy * (Y / denom) + intrinsic_.cy;

  return point2d;
}

Eigen::Vector3d EUCM::unproject(const Eigen::Vector2d & point2d) const
{
  const double fx = intrinsic_.fx;
  const double fy = intrinsic_.fy;
  const double cx = intrinsic_.cx;
  const double cy = intrinsic_.cy;
  const double alpha = distortion_.alpha;
  const double beta = distortion_.beta;
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
