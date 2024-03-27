#include "model/KB8.hpp"

namespace FCA
{
namespace model
{

KB8::KB8(const std::string & model_name, const std::string & config_path)
: Base(model_name, config_path)
{
  parse();
}

void KB8::parse() {}

void KB8::initialize(
  const Base::Params & common_params, const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
  assert(point3d_vec.size() == point2d_vec.size());

  // set fx,fy,cx,cy
  common_params_ = common_params;

  // set k1, k2, k3, k4
  Eigen::MatrixXd A(point3d_vec.size() * 2, 4);
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
    const double theta = std::atan2(r, Z);
    const double theta2 = theta * theta;
    const double theta3 = theta * theta2;
    const double theta5 = theta3 * theta2;
    const double theta7 = theta5 * theta2;
    const double theta9 = theta7 * theta2;

    A(i * 2, 0) = theta3;
    A(i * 2, 1) = theta5;
    A(i * 2, 2) = theta7;
    A(i * 2, 3) = theta9;
    A(i * 2 + 1, 0) = theta3;
    A(i * 2 + 1, 1) = theta5;
    A(i * 2 + 1, 2) = theta7;
    A(i * 2 + 1, 3) = theta9;
    b[i * 2] = (u - common_params_.cx) * (r / (common_params_.fx * X)) - theta;
    b[i * 2 + 1] = (v - common_params_.cy) * (r / (common_params_.fy * Y)) - theta;
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
  point2d.x() = common_params_.fx * d_theta * (X / r) + common_params_.cx;
  point2d.y() = common_params_.fy * d_theta * (Y / r) + common_params_.cy;

  return point2d;
}

Eigen::Vector3d KB8::unproject(const Eigen::Vector2d & point2d) const
{
  const double u = point2d.x();
  const double v = point2d.y();

  const double mx = (u - common_params_.cx) / common_params_.fx;
  const double my = (v - common_params_.cy) / common_params_.fy;

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
  }
  Eigen::Vector3d point3d;
  point3d.x() = sin(theta) * (mx / ru);
  point3d.y() = sin(theta) * (my / ru);
  point3d.z() = cos(theta);

  return point3d;
}

void KB8::optimize(
  const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
  double parameters[8] = {common_params_.fx, common_params_.fy, common_params_.cx,
                          common_params_.cy, distortion_.k1,    distortion_.k2,
                          distortion_.k3,    distortion_.k4};

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

    KB8AnalyticCostFunction * cost_function =
      new KB8AnalyticCostFunction(gt_u, gt_v, obs_x, obs_y, obs_z);
    problem.AddResidualBlock(cost_function, nullptr, parameters);
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
  distortion_.k1 = parameters[4];
  distortion_.k2 = parameters[5];
  distortion_.k3 = parameters[6];
  distortion_.k4 = parameters[7];
}

void KB8::print() const
{
  std::cout << "Final parameters: "
            << "fx=" << common_params_.fx << ", "
            << "fy=" << common_params_.fy << ", "
            << "cx=" << common_params_.cx << ", "
            << "cy=" << common_params_.cy << ", "
            << "k1=" << distortion_.k1 << ", "
            << "k2=" << distortion_.k2 << ", "
            << "k3=" << distortion_.k3 << ", "
            << "k4=" << distortion_.k4 << std::endl;
}

void KB8::save_result(const std::string & result_path) const
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
  out << YAML::Key << "coefficients" << YAML::Value;
  out << YAML::Flow << YAML::BeginSeq << distortion_.k1 << distortion_.k2 << distortion_.k3
      << distortion_.k4 << YAML::EndSeq;
  out << YAML::EndMap;

  out << YAML::EndMap;

  std::ofstream fout(result_path + "/" + model_name_ + ".yml");
  fout << out.c_str();
  fout << std::endl;
}
}  // namespace model
}  // namespace FCA
