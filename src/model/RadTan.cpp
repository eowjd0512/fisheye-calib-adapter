/*
MIT License

Copyright (c) 2024 Sangjun Lee, STRADVISION

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "model/RadTan.hpp"

namespace FCA
{
namespace model
{

RadTan::RadTan(const std::string & model_name, const std::string & config_path)
: Base(model_name, config_path)
{
}

void RadTan::parse()
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
  distortion_.p1 = config["parameter"]["p1"].as<double>();
  distortion_.p2 = config["parameter"]["p2"].as<double>();
}

void RadTan::set_sample_points(const std::vector<Eigen::Vector2d> & point2d_vec){};

void RadTan::initialize(
  const Base::Params & common_params, const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
  assert(point3d_vec.size() == point2d_vec.size());

  // set fx,fy,cx,cy
  common_params_ = common_params;

  // set k1, k2, k3
  Eigen::MatrixXd A(point3d_vec.size() * 2, 3);
  Eigen::VectorXd b(point3d_vec.size() * 2, 1);

  for (auto i = 0U; i < point3d_vec.size(); ++i) {
    const auto & point3d = point3d_vec.at(i);
    const auto & point2d = point2d_vec.at(i);
    const double X = point3d.x();
    const double Y = point3d.y();
    const double Z = point3d.z();
    const double u = point2d.x();
    const double v = point2d.y();

    const double x_prime = X / Z;
    const double y_prime = Y / Z;

    const double r2 = x_prime * x_prime + y_prime * y_prime;
    const double r4 = r2 * r2;
    const double r6 = r4 * r2;

    A(i * 2, 0) = r2;
    A(i * 2, 1) = r4;
    A(i * 2, 2) = r6;
    A(i * 2 + 1, 0) = r2;
    A(i * 2 + 1, 1) = r4;
    A(i * 2 + 1, 2) = r6;
    b[i * 2] = (u - common_params_.cx) / (common_params_.fx * x_prime) - 1.0;
    b[i * 2 + 1] = (v - common_params_.cy) / (common_params_.fy * y_prime) - 1.0;
  }

  const Eigen::VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
  distortion_.k1 = x[0];
  distortion_.k2 = x[1];
  distortion_.k3 = x[2];
  distortion_.p1 = 0.0;
  distortion_.p2 = 0.0;
}

Eigen::Vector2d RadTan::project(const Eigen::Vector3d & point3d, bool condition) const
{
  const double X = point3d.x();
  const double Y = point3d.y();
  const double Z = point3d.z();
  const double fx = common_params_.fx;
  const double fy = common_params_.fy;
  const double cx = common_params_.cx;
  const double cy = common_params_.cy;

  const double k1 = distortion_.k1;
  const double k2 = distortion_.k2;
  const double k3 = distortion_.k3;
  const double p1 = distortion_.p1;
  const double p2 = distortion_.p2;

  const double x_prime = X / Z;
  const double y_prime = Y / Z;

  const double r2 = x_prime * x_prime + y_prime * y_prime;
  const double r4 = r2 * r2;
  const double r6 = r4 * r2;

  const double x_distorted_estimate = x_prime * (1.0 + k1 * r2 + k2 * r4 + k3 * r6) +
                                      2.0 * p1 * x_prime * y_prime +
                                      p2 * (r2 + 2.0 * x_prime * x_prime);
  const double y_distorted_estimate = y_prime * (1.0 + k1 * r2 + k2 * r4 + k3 * r6) +
                                      p1 * (r2 + 2.0 * y_prime * y_prime) +
                                      2.0 * p2 * x_prime * y_prime;

  Eigen::Vector2d point2d;
  point2d.x() = fx * x_distorted_estimate + cx;
  point2d.y() = fy * y_distorted_estimate + cy;

  return point2d;
}

Eigen::Vector3d RadTan::unproject(const Eigen::Vector2d & point2d, bool condition) const
{
  const double u = point2d.x();
  const double v = point2d.y();
  const double fx = common_params_.fx;
  const double fy = common_params_.fy;
  const double cx = common_params_.cx;
  const double cy = common_params_.cy;

  const double k1 = distortion_.k1;
  const double k2 = distortion_.k2;
  const double k3 = distortion_.k3;
  const double p1 = distortion_.p1;
  const double p2 = distortion_.p2;

  const double x_distorted = (u - cx) / fx;
  const double y_distorted = (v - cy) / fy;

  Eigen::Vector2d point(x_distorted, y_distorted);

  constexpr auto EPS = 1e-6;

  while (true) {
    const double x = point[0];
    const double y = point[1];
    const double r2 = x * x + y * y;
    const double r4 = r2 * r2;
    const double r6 = r4 * r2;

    const double d = (1.0 + k1 * r2 + k2 * r4 + k3 * r6);

    const double x_distorted_estimate = x * d + 2.0 * p1 * x * y + p2 * (r2 + 2.0 * x * x);
    const double y_distorted_estimate = y * d + p1 * (r2 + 2.0 * y * y) + 2.0 * p2 * x * y;

    Eigen::Vector2d error(x_distorted_estimate - x_distorted, y_distorted_estimate - y_distorted);

    if (error.norm() < EPS) {
      break;
    }

    Eigen::Matrix2d J;
    const double term1 = (k1 + 2.0 * r2 * k2 + 3.0 * r4 * k3);
    J(0, 0) = d + 2.0 * x * term1 + 2.0 * p1 * y + 6.0 * p2 * x;
    J(0, 1) = 2.0 * x * y * term1 + 2.0 * p1 * x + 2.0 * p2 * y;
    J(1, 0) = 2.0 * y * x * term1 + 2.0 * p2 * y + 2.0 * p1 * x;
    J(1, 1) = d + 2.0 * y * term1 + 2.0 * p2 * x + 6.0 * p1 * y;

    Eigen::Vector2d delta = J.inverse() * error;
    point -= delta;

    if (delta.norm() < EPS) {
      break;
    }
  }

  Eigen::Vector3d point3d;
  point3d.x() = point.x();
  point3d.y() = point.y();
  point3d.z() = 1.0;
  point3d.normalized();
  return point3d;
}

bool RadTan::check_proj_condition(double z) { return z > 0.0; }

void RadTan::optimize(
  const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec, bool display_optimization_progress)
{
  double parameters[9] = {common_params_.fx, common_params_.fy, common_params_.cx,
                          common_params_.cy, distortion_.k1,    distortion_.k2,
                          distortion_.k3,    distortion_.p1,    distortion_.p2};

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

    RadTanAnalyticCostFunction * cost_function =
      new RadTanAnalyticCostFunction(gt_u, gt_v, obs_x, obs_y, obs_z);
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
  distortion_.p1 = parameters[7];
  distortion_.p2 = parameters[8];
}

void RadTan::print() const
{
  std::cout << model_name_ << " parameters: "
            << "fx=" << common_params_.fx << ", "
            << "fy=" << common_params_.fy << ", "
            << "cx=" << common_params_.cx << ", "
            << "cy=" << common_params_.cy << ", "
            << "k1=" << distortion_.k1 << ", "
            << "k2=" << distortion_.k2 << ", "
            << "k3=" << distortion_.k3 << ", "
            << "p1=" << distortion_.p1 << ", "
            << "p2=" << distortion_.p2 << std::endl;
}

void RadTan::save_result(const std::string & result_path) const
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
  out << YAML::Key << "p1" << YAML::Value << distortion_.p1;
  out << YAML::Key << "p2" << YAML::Value << distortion_.p2;
  out << YAML::EndMap;

  out << YAML::EndMap;

  std::ofstream fout(result_path + "/" + model_name_ + ".yml");
  fout << out.c_str();
  fout << std::endl;
}

void RadTan::evaluate(const model::Base * const gt)
{
  const RadTan * gt_model = dynamic_cast<const RadTan *>(gt);

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
  const double diff_p1 = est_distortion_params.p1 - gt_distortion_params.p1;
  const double diff_p2 = est_distortion_params.p2 - gt_distortion_params.p2;

  const double params_diff_norm = std::sqrt(
    diff_fx * diff_fx + diff_fy * diff_fy + diff_cx * diff_cx + diff_cy * diff_cy +
    diff_k1 * diff_k1 + diff_k2 * diff_k2 + diff_k3 * diff_k3 + diff_p1 * diff_p1 +
    diff_p2 * diff_p2);

  std::cout << "parameter error: " << params_diff_norm << std::endl;
}
}  // namespace model
}  // namespace FCA
