#include "model/OcamLib.hpp"

#include <limits>

namespace FCA
{
namespace model
{

OcamLib::OcamLib(const std::string & model_name, const std::string & config_path)
: Base(model_name, config_path)
{
  parse();
}

void OcamLib::parse()
{
  // Load the YAML file
  const YAML::Node config = YAML::LoadFile(config_path_ + "/" + model_name_ + ".yml");

  // Read coefficients
  distortion_.proj_coeffs = config["coefficients"]["projection"].as<std::vector<double>>();
  distortion_.unproj_coeffs = config["coefficients"]["unprojection"].as<std::vector<double>>();

  // Read parameters
  const std::vector<double> affine_parameters =
    config["parameter"]["affine"].as<std::vector<double>>();
  distortion_.c = affine_parameters[0];
  distortion_.d = affine_parameters[1];
  distortion_.e = affine_parameters[2];
  intrinsic_.fx = distortion_.proj_coeffs.at(0);  // approximation
  intrinsic_.fy = distortion_.proj_coeffs.at(0);  // approximation
  intrinsic_.cx = config["parameter"]["cx"].as<double>();
  intrinsic_.cy = config["parameter"]["cy"].as<double>();
}

void OcamLib::initialize(
  const Base::Params & intrinsic, const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
  assert(point3d_vec.size() == point2d_vec.size());

  // set fx,fy,cx,cy
  intrinsic_ = intrinsic;

  // set c, d, e
  distortion_.c = 1.0;
  distortion_.d = 0.0;
  distortion_.e = 0.0;

  // set polynomial coefficients for unprojection
  Eigen::MatrixXd A(point3d_vec.size(), 5);
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
    const double r2 = r * r;
    const double r3 = r * r2;
    const double r4 = r2 * r2;

    A(i, 0) = 1;
    A(i, 1) = r;
    A(i, 2) = r2;
    A(i, 3) = r3;
    A(i, 4) = r4;
    b[i] = (u - intrinsic_.cx) * (Z / X);
  }

  const Eigen::VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
  distortion_.unproj_coeffs = std::vector<double>(x.data(), x.data() + x.size());
}

Eigen::Vector2d OcamLib::project(const Eigen::Vector3d & point3d) const
{
  const double X = point3d.x();
  const double Y = point3d.y();
  const double Z = point3d.z();
  const double r = std::sqrt((X * X) + (Y * Y));
  const double theta = atan(Z / r);

  double rho = distortion_.unproj_coeffs.at(0);
  double theta_i = 1.0;

  for (auto i = 1U; i < distortion_.unproj_coeffs.size(); ++i) {
    theta_i *= theta;
    rho += theta_i * distortion_.unproj_coeffs.at(i);
  }

  const double u = (X / r) * rho;
  const double v = (Y / r) * rho;

  Eigen::Vector2d point2d;
  point2d.x() = (u * distortion_.c) + (v * distortion_.d) + intrinsic_.cx;
  point2d.y() = (u * distortion_.e) + v + intrinsic_.cy;

  return point2d;
}

void OcamLib::estimate_projection_coefficients(
  const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
  assert(point3d_vec.size() == point2d_vec.size());

  constexpr auto MAX_POLY_NUM = 20;

  std::map<uint32_t, std::vector<double>> coefficient_candidates;
  double max_error = std::numeric_limits<double>::max();
  int32_t best_poly_num = -1;
  for (auto n = 0; n <= MAX_POLY_NUM; ++n) {
    // set polynomial coefficients for projection
    Eigen::MatrixXd A(point3d_vec.size() * 2, n + 1);
    Eigen::VectorXd b(point3d_vec.size() * 2, 1);

    for (auto i = 0U; i < point3d_vec.size(); ++i) {
      const auto & point3d = point3d_vec.at(i);
      const auto & point2d = point2d_vec.at(i);
      const double X = point3d.x();
      const double Y = point3d.y();
      const double Z = point3d.z();
      const double u = point2d.x();
      const double v = point2d.y();
      const double r = std::sqrt((X * X) + (Y * Y));
      const double theta = atan(Z / r);

      const double term1 = (distortion_.c * (X / r)) + (distortion_.d * (Y / r));
      const double term2 = (distortion_.e * (X / r)) + (Y / r);

      double theta_i = 1.0;
      A(2 * i, 0) = term1;
      A((2 * i) + 1, 0) = term2;
      for (auto j = 1; j <= n; ++j) {
        theta_i *= theta;
        A(2 * i, j) = term1 * theta_i;
        A((2 * i) + 1, j) = term2 * theta_i;
      }
      b[2 * i] = u - intrinsic_.cx;
      b[(2 * i) + 1] = v - intrinsic_.cy;
    }
    const Eigen::VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
    distortion_.proj_coeffs.clear();
    distortion_.proj_coeffs = std::vector<double>(x.data(), x.data() + x.size());

    const double error = evaluate(point3d_vec, point2d_vec);
    if (error < max_error) {
      max_error = error;
      coefficient_candidates.insert(std::make_pair(n, distortion_.proj_coeffs));
      best_poly_num = n;
    }
  }
  distortion_.proj_coeffs.clear();
  distortion_.proj_coeffs = coefficient_candidates.at(best_poly_num);
}

Eigen::Vector3d OcamLib::unproject(const Eigen::Vector2d & point2d) const
{
  const double u = point2d.x();
  const double v = point2d.y();
  // 1/det(A), where A = [c,d;e,1]
  const double inv_det = 1.0 / (distortion_.c - (distortion_.d * distortion_.e));

  const double px = inv_det * ((u - intrinsic_.cx) - distortion_.d * (v - intrinsic_.cy));
  const double py =
    inv_det * (-distortion_.e * (u - intrinsic_.cx) + distortion_.c * (v - intrinsic_.cy));

  //distance [pixels] of  the point from the image center
  const double r = std::sqrt((px * px) + (py * py));

  double pz = distortion_.proj_coeffs.at(0);
  double r_i = 1.0;
  for (auto i = 1U; i < distortion_.proj_coeffs.size(); ++i) {
    r_i *= r;
    pz += r_i * distortion_.proj_coeffs.at(i);
  }

  const double norm = std::sqrt((px * px) + (py * py) + (pz * pz));

  Eigen::Vector3d point3d;
  point3d.x() = px / norm;
  point3d.y() = py / norm;
  point3d.z() = pz / norm;

  return point3d;
}

Eigen::MatrixXd OcamLib::calculate_jacobian() const {}

void OcamLib::optimize() {}

void OcamLib::print() const {}

double OcamLib::evaluate(
  const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
  double distance_sum = 0.0;
  for (auto i = 0U; i < point3d_vec.size(); ++i) {
    const auto & point3d = point3d_vec.at(i);
    const auto & point2d = point2d_vec.at(i);

    const Eigen::Vector2d estimated_point2d = project(point3d);
    distance_sum += (point2d - estimated_point2d).norm();
  }
  return distance_sum / point3d_vec.size();
}

}  // namespace model
}  // namespace FCA
