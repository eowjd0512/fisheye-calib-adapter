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
  const Base::Params & common_params, const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
  assert(point3d_vec.size() == point2d_vec.size());

  // set fx,fy,cx,cy
  common_params_ = common_params;

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
  const double pc = u - common_params_.cx;
  distortion_.alpha = ((common_params_.fx * X) - (pc * Z)) / (pc * (d - Z));
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
  point2d.x() = common_params_.fx * (X / denom) + common_params_.cx;
  point2d.y() = common_params_.fy * (Y / denom) + common_params_.cy;

  return point2d;
}

Eigen::Vector3d EUCM::unproject(const Eigen::Vector2d & point2d) const
{
  const double fx = common_params_.fx;
  const double fy = common_params_.fy;
  const double cx = common_params_.cx;
  const double cy = common_params_.cy;
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

void EUCM::optimize()
{
  double parameters[6] = {common_params_.fx, common_params_.fy, common_params_.cx,
                          common_params_.cy, distortion_.alpha, distortion_.beta};

  double gt_u = 0.5;
  double gt_v = 0.5;
  double obs_x = 0.5;
  double obs_y = 0.5;
  double obs_z = 0.5;
  ceres::Problem problem;
  EUCMAnalyticCostFunction * cost_function =
    new EUCMAnalyticCostFunction(gt_u, gt_v, obs_x, obs_y, obs_z);
  problem.AddResidualBlock(cost_function, nullptr, parameters);

  //// Auto diff
  // problem.AddResidualBlock(
  //     new ceres::AutoDiffCostFunction<EUCMAutoDiffCostFunctor, 2, 6>(
  //         new EUCMAutoDiffCostFunctor(observed_x, observed_y)),
  //     nullptr, parameters);

  //// set parameters range
  // problem.SetParameterLowerBound(parameters, 0, 0.1); // fx > 0.1
  //   problem.SetParameterLowerBound(parameters, 1, 0.1); // fy > 0.1
  //   problem.SetParameterLowerBound(parameters, 4, 0.0); // alpha >= 0
  //   problem.SetParameterUpperBound(parameters, 4, 1.0); // alpha <= 1
  //   problem.SetParameterLowerBound(parameters, 5, 0.0); // beta >= 0
  //   problem.SetParameterUpperBound(parameters, 5, 1.0); // beta <= 1

  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = true;

  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  std::cout << summary.FullReport() << "\n";
  std::cout << "Final parameters: "
            << "fx=" << parameters[0] << ", "
            << "fy=" << parameters[1] << ", "
            << "cx=" << parameters[2] << ", "
            << "cy=" << parameters[3] << ", "
            << "alpha=" << parameters[4] << ", "
            << "beta=" << parameters[5] << std::endl;
}

void EUCM::print() const {}

}  // namespace model
}  // namespace FCA
