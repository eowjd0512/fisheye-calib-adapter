#include "model/KB.hpp"

namespace FCA
{
namespace model
{

KB::KB(const std::string & model_name, const std::string & config_path)
: Base(model_name, config_path)
{
}

void KB::parse()
{
  // Load the YAML file
  YAML::Node config;
  try {
    config = YAML::LoadFile(config_path_ + "/" + model_name_ + ".yml");
  } catch (const YAML::BadFile & e) {
    std::cerr << "Error loading YAML file: " << e.what() << std::endl;
    exit(-1);
  }

  common_params_.width = config["image"]["width"].as<int32_t>();
  common_params_.height = config["image"]["height"].as<int32_t>();

  // Read parameters
  common_params_.cx = config["parameter"]["cx"].as<double>();
  common_params_.cy = config["parameter"]["cy"].as<double>();
  common_params_.fx = config["parameter"]["fx"].as<double>();
  common_params_.fy = config["parameter"]["fy"].as<double>();
  distortion_.k1 = config["parameter"]["k1"].as<double>();
  distortion_.k2 = config["parameter"]["k2"].as<double>();
  distortion_.k3 = config["parameter"]["k3"].as<double>();
  distortion_.k4 = config["parameter"]["k4"].as<double>();
}

void KB::set_sample_points(const std::vector<Eigen::Vector2d> & point2d_vec){};

void KB::initialize(
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

Eigen::Vector2d KB::project(const Eigen::Vector3d & point3d, bool condition) const
{
  const double X = point3d.x();
  const double Y = point3d.y();
  const double Z = point3d.z();

  if (condition) {
    if (!check_proj_condition(Z)) {
      return {-1., -1.};
    }
  }

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

Eigen::Vector3d KB::unproject(const Eigen::Vector2d & point2d, bool condition) const
{
  const double u = point2d.x();
  const double v = point2d.y();

  const double mx = (u - common_params_.cx) / common_params_.fx;
  const double my = (v - common_params_.cy) / common_params_.fy;

  double ru = std::sqrt(mx * mx + my * my);
  ru = std::min(std::max(-M_PI / 2.0, ru), M_PI / 2.0);
  double theta = ru;
  constexpr double PRECISION = 1e-3;

  bool converged = true;
  constexpr auto MAX_ITERATION = 10;
  auto i = 0;
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
      } else if (++i > MAX_ITERATION) {
        converged = false;
        break;
      }
    }
  } else {
    converged = false;
  }

  if (!converged) {
    return {-1., -1., -1.};
  }

  Eigen::Vector3d point3d;
  point3d.x() = sin(theta) * (mx / ru);
  point3d.y() = sin(theta) * (my / ru);
  point3d.z() = cos(theta);

  return point3d;
}

bool KB::check_proj_condition(double z) { return z > 0.0; }

void KB::optimize(
  const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec, bool display_optimization_progress)
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

    if (!check_proj_condition(obs_z)) {
      continue;
    }
    KBAnalyticCostFunction * cost_function =
      new KBAnalyticCostFunction(gt_u, gt_v, obs_x, obs_y, obs_z);
    problem.AddResidualBlock(cost_function, nullptr, parameters);
  }
  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = false;
  if (display_optimization_progress) {
    options.minimizer_progress_to_stdout = true;
  }

  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  if (display_optimization_progress) {
    std::cout << summary.FullReport() << std::endl;
  }

  common_params_.fx = parameters[0];
  common_params_.fy = parameters[1];
  common_params_.cx = parameters[2];
  common_params_.cy = parameters[3];
  distortion_.k1 = parameters[4];
  distortion_.k2 = parameters[5];
  distortion_.k3 = parameters[6];
  distortion_.k4 = parameters[7];
}

void KB::print() const
{
  std::cout << model_name_ << " parameters: "
            << "fx=" << common_params_.fx << ", "
            << "fy=" << common_params_.fy << ", "
            << "cx=" << common_params_.cx << ", "
            << "cy=" << common_params_.cy << ", "
            << "k1=" << distortion_.k1 << ", "
            << "k2=" << distortion_.k2 << ", "
            << "k3=" << distortion_.k3 << ", "
            << "k4=" << distortion_.k4 << std::endl;
}

void KB::save_result(const std::string & result_path) const
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
  out << YAML::Key << "k1" << YAML::Value << distortion_.k1;
  out << YAML::Key << "k2" << YAML::Value << distortion_.k2;
  out << YAML::Key << "k3" << YAML::Value << distortion_.k3;
  out << YAML::Key << "k4" << YAML::Value << distortion_.k4;
  out << YAML::EndMap;

  out << YAML::EndMap;

  std::ofstream fout(result_path + "/" + model_name_ + ".yml");
  fout << out.c_str();
  fout << std::endl;
}

void KB::evaluate(const model::Base * const gt)
{
  const KB * gt_model = dynamic_cast<const KB *>(gt);

  const auto & est_pinhole_params = this->common_params_;
  const auto & gt_pinhole_params = gt_model->get_common_params();
  const auto & est_distortion_params = this->distortion_;
  const auto & gt_distortion_params = gt_model->get_distortion_params();

  const double diff_fx = est_pinhole_params.fx - gt_pinhole_params.fx;
  const double diff_fy = est_pinhole_params.fy - gt_pinhole_params.fy;
  const double diff_cx = est_pinhole_params.cx - gt_pinhole_params.cx;
  const double diff_cy = est_pinhole_params.cy - gt_pinhole_params.cy;
  const double diff_k1 = est_distortion_params.k1 - gt_distortion_params.k1;
  const double diff_k2 = est_distortion_params.k2 - gt_distortion_params.k2;
  const double diff_k3 = est_distortion_params.k3 - gt_distortion_params.k3;
  const double diff_k4 = est_distortion_params.k4 - gt_distortion_params.k4;

  const double params_diff_norm = std::sqrt(
    diff_fx * diff_fx + diff_fy * diff_fy + diff_cx * diff_cx + diff_cy * diff_cy +
    diff_k1 * diff_k1 + diff_k2 * diff_k2 + diff_k3 * diff_k3 + diff_k4 * diff_k4);

  std::cout << "parameter error: " << params_diff_norm << std::endl;
}
}  // namespace model
}  // namespace FCA
