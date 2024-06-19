#include "model/double_sphere.hpp"

namespace FCA
{
namespace model
{

DoubleSphere::DoubleSphere(const std::string & model_name, const std::string & config_path)
: Base(model_name, config_path)
{
  parse();
}

void DoubleSphere::parse()
{
  // Load the YAML file
  const YAML::Node config = YAML::LoadFile(config_path_ + "/" + model_name_ + ".yml");

  common_params_.width = config["image"]["width"].as<int32_t>();
  common_params_.height = config["image"]["height"].as<int32_t>();

  // Read parameters
  common_params_.cx = config["parameter"]["cx"].as<double>();
  common_params_.cy = config["parameter"]["cy"].as<double>();
  common_params_.fx = config["parameter"]["fx"].as<double>();
  common_params_.fy = config["parameter"]["fy"].as<double>();
  distortion_.alpha = config["parameter"]["alpha"].as<double>();
  distortion_.xi = config["parameter"]["xi"].as<double>();
}

void DoubleSphere::initialize(
  const Base::Params & common_params, const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
  assert(point3d_vec.size() == point2d_vec.size());

  // set fx,fy,cx,cy
  common_params_ = common_params;

  // set xi
  distortion_.xi = 0.0;

  // set alpha
  Eigen::MatrixXd A(point3d_vec.size() * 2, 1);
  Eigen::VectorXd b(point3d_vec.size() * 2, 1);

  for (auto i = 0U; i < point3d_vec.size(); ++i) {
    const auto & point3d = point3d_vec.at(i);
    const auto & point2d = point2d_vec.at(i);
    const double X = point3d.x();
    const double Y = point3d.y();
    const double Z = point3d.z();
    const double u = point2d.x();
    const double v = point2d.y();

    const double d = std::sqrt((X * X) + (Y * Y) + (Z * Z));
    const double u_cx = u - common_params_.cx;
    const double v_cy = v - common_params_.cy;
    // distortion_.alpha = ((common_params_.fx * X) - (u_cx * Z)) / (u_cx * (d - Z));

    A(i * 2, 0) = u_cx * (d - Z);
    A(i * 2 + 1, 0) = v_cy * (d - Z);

    b[i * 2] = (common_params_.fx * X) - (u_cx * Z);
    b[i * 2 + 1] = (common_params_.fy * Y) - (v_cy * Z);
  }

  const Eigen::VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
  distortion_.alpha = x[0];

  std::cout << "initialized alpha: " << distortion_.alpha << std::endl;
}

Eigen::Vector2d DoubleSphere::project(const Eigen::Vector3d & point3d) const
{
  const double X = point3d.x();
  const double Y = point3d.y();
  const double Z = point3d.z();

  const double r_squared = (X * X) + (Y * Y);
  const double d1 = std::sqrt(r_squared + (Z * Z));
  const double gamma = distortion_.xi * d1 + Z;
  const double d2 = std::sqrt(r_squared + gamma * gamma);

  const double denom = distortion_.alpha * d2 + (1.0 - distortion_.alpha) * gamma;

  constexpr double PRECISION = 1e-3;
  if ((denom < PRECISION) || !check_proj_condition(Z, d1, distortion_.xi, distortion_.alpha)) {
    return {-1, -1};
  }

  Eigen::Vector2d point2d;
  point2d.x() = common_params_.fx * (X / denom) + common_params_.cx;
  point2d.y() = common_params_.fy * (Y / denom) + common_params_.cy;

  return point2d;
}

Eigen::Vector3d DoubleSphere::unproject(const Eigen::Vector2d & point2d) const
{
  const double fx = common_params_.fx;
  const double fy = common_params_.fy;
  const double cx = common_params_.cx;
  const double cy = common_params_.cy;
  const double alpha = distortion_.alpha;
  const double xi = distortion_.xi;

  const double u = point2d.x();
  const double v = point2d.y();
  const double gamma = 1.0 - alpha;
  const double mx = (u - cx) / fx * gamma;
  const double my = (v - cy) / fy * gamma;
  const double r_squared = (mx * mx) + (my * my);
  const double mz = (1.0 - alpha * alpha * r_squared) /
                    (alpha * std::sqrt(1.0 - (2.0 * alpha - 1.0) * r_squared) + gamma);
  const double mz_squared = mz * mz;

  const double num = mz * xi + std::sqrt(mz_squared + (1.0 - xi * xi) * r_squared);
  const double denom = mz_squared + r_squared;

  constexpr double PRECISION = 1e-3;
  if ((denom < PRECISION) || !check_unproj_condition(r_squared, alpha)) {
    return {-1, -1, -1};
  }

  const double coeff = num / denom;

  Eigen::Vector3d point3d;
  point3d = coeff * Eigen::Vector3d(mx, my, mz) - Eigen::Vector3d(0., 0., xi);

  return point3d;
}

bool DoubleSphere::check_proj_condition(double z, double d1, double xi, double alpha)
{
  double w1 = (1.0 - alpha) / alpha;
  if (alpha <= 0.5) {
    w1 = alpha / (1.0 - alpha);
  }
  const double w2 = (w1 + xi) / std::sqrt(2 * w1 * xi + xi * xi + 1.0);
  return z > -w2 * d1;
}

bool DoubleSphere::check_unproj_condition(double r_squared, double alpha)
{
  bool condition = true;
  if (alpha > 0.5) {
    if (r_squared > 1.0 / (2 * alpha - 1.0)) {
      condition = false;
    }
  }
  return condition;
}

void DoubleSphere::optimize(
  const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
  double parameters[6] = {common_params_.fx, common_params_.fy, common_params_.cx,
                          common_params_.cy, distortion_.alpha, distortion_.xi};

  ceres::Problem problem;

  const auto num_pairs = point3d_vec.size();
  for (auto i = 0U; i < num_pairs; ++i) {
    const auto & point2d = point2d_vec.at(i);
    const auto & point3d = point3d_vec.at(i);
    const double gt_u = point2d.x();
    const double gt_v = point2d.y();
    const double obs_x = point3d.x();
    const double obs_y = point3d.y();
    const double obs_z = point3d.z();

    DSAnalyticCostFunction * cost_function =
      new DSAnalyticCostFunction(gt_u, gt_v, obs_x, obs_y, obs_z);
    problem.AddResidualBlock(cost_function, nullptr, parameters);

    // set parameters range
    problem.SetParameterLowerBound(parameters, 4, 0.0001);  // alpha >= 0
    problem.SetParameterUpperBound(parameters, 4, 0.999);  // alpha <= 1
  }
  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = true;

  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  std::cout << summary.FullReport() << std::endl;

  common_params_.fx = parameters[0];
  common_params_.fy = parameters[1];
  common_params_.cx = parameters[2];
  common_params_.cy = parameters[3];
  distortion_.alpha = parameters[4];
  distortion_.xi = parameters[5];
}

void DoubleSphere::print() const
{
  std::cout << "Final parameters: "
            << "fx=" << common_params_.fx << ", "
            << "fy=" << common_params_.fy << ", "
            << "cx=" << common_params_.cx << ", "
            << "cy=" << common_params_.cy << ", "
            << "alpha=" << distortion_.alpha << ", "
            << "xi=" << distortion_.xi << std::endl;
}

void DoubleSphere::save_result(const std::string & result_path) const
{
  YAML::Emitter out;

  out << YAML::BeginMap;
  out << YAML::Key << "image" << YAML::Value;
  out << YAML::BeginMap;
  out << YAML::Key << "width" << YAML::Value << common_params_.width;
  out << YAML::Key << "height" << YAML::Value << common_params_.height;
  out << YAML::EndMap;

  out << YAML::Key << "parameter" << YAML::Value;
  out << YAML::BeginMap;
  out << YAML::Key << "fx" << YAML::Value << common_params_.fx;
  out << YAML::Key << "fy" << YAML::Value << common_params_.fy;
  out << YAML::Key << "cx" << YAML::Value << common_params_.cx;
  out << YAML::Key << "cy" << YAML::Value << common_params_.cy;
  out << YAML::Key << "alpha" << YAML::Value << distortion_.alpha;
  out << YAML::Key << "xi" << YAML::Value << distortion_.xi;
  out << YAML::EndMap;

  out << YAML::EndMap;

  std::ofstream fout(result_path + "/" + model_name_ + ".yml");
  fout << out.c_str();
  fout << std::endl;
}
}  // namespace model
}  // namespace FCA
