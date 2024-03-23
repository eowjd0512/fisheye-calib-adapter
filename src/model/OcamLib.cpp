#include "model/OcamLib.hpp"

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
  params.proj_coeffs = config["coefficients"]["projection"].as<std::vector<double>>();
  params.unproj_coeffs = config["coefficients"]["unprojection"].as<std::vector<double>>();

  // Read parameters
  const std::vector<double> affine_parameters =
    config["parameter"]["affine"].as<std::vector<double>>();
  params.c = affine_parameters[0];
  params.d = affine_parameters[1];
  params.e = affine_parameters[2];
  params.cx = config["parameter"]["cx"].as<double>();
  params.cy = config["parameter"]["cy"].as<double>();
}

void OcamLib::initialize() {}

Eigen::Vector2d OcamLib::project(const Eigen::Vector3d & point3d) const
{
  const double X = point3d.x();
  const double Y = point3d.y();
  const double Z = point3d.z();
  const double r = std::sqrt(X * X + Y * Y);
  const double theta = atan(Z / r);

  double rho = params.unproj_coeffs[0];
  double theta_i = 1.0;

  for (auto i = 1U; i < params.unproj_coeffs.size(); ++i) {
    theta_i *= theta;
    rho += theta_i * params.unproj_coeffs[i];
  }

  const double u = X / r * rho;
  const double v = Y / r * rho;

  Eigen::Vector2d point2d;
  point2d.x() = u * params.c + v * params.d + params.cx;
  point2d.y() = u * params.e + v + params.cy;

  return point2d;
}

Eigen::Vector3d OcamLib::unproject(const Eigen::Vector2d & point2d) const
{
  const double u = point2d.x();
  const double v = point2d.y();
  // 1/det(A), where A = [c,d;e,1]
  const double inv_det = 1.0 / (params.c - params.d * params.e);

  const double px = inv_det * ((u - params.cx) - params.d * (v - params.cy));
  const double py = inv_det * (-params.e * (u - params.cx) + params.c * (v - params.cy));

  //distance [pixels] of  the point from the image center
  const double r = std::sqrt(px * px + py * py);

  double pz = params.proj_coeffs[0];
  double r_i = 1.0;
  for (auto i = 1U; i < params.proj_coeffs.size(); ++i) {
    r_i *= r;
    pz += r_i * params.proj_coeffs[i];
  }

  const double norm = std::sqrt(px * px + py * py + pz * pz);

  Eigen::Vector3d point3d;
  point3d.x() = px / norm;
  point3d.y() = py / norm;
  point3d.z() = pz / norm;

  return point3d;
}

Eigen::MatrixXd OcamLib::calculate_jacobian() const {}

void OcamLib::optimize() {}

void OcamLib::print() const {}

}  // namespace model
}  // namespace FCA
