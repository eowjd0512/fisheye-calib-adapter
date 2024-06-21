#include "model/OcamLib.hpp"

#include <limits>

namespace FCA
{
namespace model
{

OcamLib::OcamLib(const std::string & model_name, const std::string & config_path)
: Base(model_name, config_path)
{
}

void OcamLib::parse()
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

  // Read coefficients
  distortion_.proj_coeffs =
    config["parameter"]["coefficients"]["projection"].as<std::vector<double>>();
  distortion_.unproj_coeffs =
    config["parameter"]["coefficients"]["unprojection"].as<std::vector<double>>();

  // Read parameters
  distortion_.c = config["parameter"]["c"].as<double>();
  distortion_.d = config["parameter"]["d"].as<double>();
  distortion_.e = config["parameter"]["e"].as<double>();
  common_params_.fx = common_params_.fy = common_params_.width / 2.0;  // approximation
  common_params_.cx = config["parameter"]["cx"].as<double>();
  common_params_.cy = config["parameter"]["cy"].as<double>();
}

void OcamLib::set_sample_points(const std::vector<Eigen::Vector2d> & point2d_vec)
{
  point2d_vec_ = point2d_vec;
};

void OcamLib::initialize(
  const Base::Params & common_params, const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
  assert(point3d_vec.size() == point2d_vec.size());

  // set fx,fy,cx,cy
  common_params_ = common_params;

  // set c, d, e
  distortion_.c = 1.0;
  distortion_.d = 0.0;
  distortion_.e = 0.0;

  // set polynomial coefficients for unprojection
  Eigen::MatrixXd A(point3d_vec.size() * 2, 5);
  Eigen::VectorXd b(point3d_vec.size() * 2, 1);

  for (auto i = 0U; i < point3d_vec.size(); ++i) {
    const auto & point3d = point3d_vec.at(i);
    const auto & point2d = point2d_vec.at(i);
    const double X = point3d.x();
    const double Y = point3d.y();
    const double Z = point3d.z();
    const double u = point2d.x();
    const double v = point2d.y();
    const double mx = u - common_params_.cx;
    const double my = v - common_params_.cy;

    const double r = std::sqrt((mx * mx) + (my * my));
    const double r2 = r * r;
    const double r3 = r * r2;
    const double r4 = r2 * r2;

    A(2 * i, 0) = 1;
    A(2 * i, 1) = r;
    A(2 * i, 2) = r2;
    A(2 * i, 3) = r3;
    A(2 * i, 4) = r4;
    A(2 * i + 1, 0) = 1;
    A(2 * i + 1, 1) = r;
    A(2 * i + 1, 2) = r2;
    A(2 * i + 1, 3) = r3;
    A(2 * i + 1, 4) = r4;
    b[2 * i] = -(u - common_params_.cx) * (Z / X);
    b[2 * i + 1] = -(v - common_params_.cy) * (Z / Y);
  }

  const Eigen::VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
  distortion_.unproj_coeffs = std::vector<double>(x.data(), x.data() + x.size());
}

Eigen::Vector2d OcamLib::project(const Eigen::Vector3d & point3d, bool condition) const
{
  const double X = point3d.x();
  const double Y = point3d.y();
  const double Z = point3d.z();
  const double r = std::sqrt((X * X) + (Y * Y));
  const double theta = atan(-Z / r);

  double rho = distortion_.proj_coeffs.at(0);
  double theta_i = 1.0;

  for (auto i = 1U; i < distortion_.proj_coeffs.size(); ++i) {
    theta_i *= theta;
    rho += theta_i * distortion_.proj_coeffs.at(i);
  }

  const double u = (X / r) * rho;
  const double v = (Y / r) * rho;

  Eigen::Vector2d point2d;
  point2d.x() = (u * distortion_.c) + (v * distortion_.d) + common_params_.cx;
  point2d.y() = (u * distortion_.e) + v + common_params_.cy;

  return point2d;
}

Eigen::Vector3d OcamLib::unproject(const Eigen::Vector2d & point2d, bool condition) const
{
  const double u = point2d.x();
  const double v = point2d.y();
  // 1/det(A), where A = [c,d;e,1]
  const double inv_det = 1.0 / (distortion_.c - (distortion_.d * distortion_.e));

  const double mx = inv_det * ((u - common_params_.cx) - distortion_.d * (v - common_params_.cy));
  const double my =
    inv_det * (-distortion_.e * (u - common_params_.cx) + distortion_.c * (v - common_params_.cy));

  //distance [pixels] of  the point from the image center
  const double r = std::sqrt((mx * mx) + (my * my));

  double mz = distortion_.unproj_coeffs.at(0);
  double r_i = 1.0;
  for (auto i = 1U; i < distortion_.unproj_coeffs.size(); ++i) {
    r_i *= r;
    mz += r_i * distortion_.unproj_coeffs.at(i);
  }

  const double norm = std::sqrt((mx * mx) + (my * my) + (mz * mz));

  Eigen::Vector3d point3d;
  point3d.x() = mx / norm;
  point3d.y() = my / norm;
  point3d.z() = -mz / norm;

  return point3d;
}

bool OcamLib::check_proj_condition(double z) { return z > 0.0; }

void OcamLib::optimize(
  const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec, bool display_optimization_progress)
{
  double parameters[7] = {
    common_params_.cx,
    common_params_.cy,
    distortion_.unproj_coeffs[0],
    distortion_.unproj_coeffs[1],
    distortion_.unproj_coeffs[2],
    distortion_.unproj_coeffs[3],
    distortion_.unproj_coeffs[4]};

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

    OcamLibAnalyticCostFunction * cost_function =
      new OcamLibAnalyticCostFunction(gt_u, gt_v, obs_x, obs_y, obs_z);
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

  common_params_.cx = parameters[0];
  common_params_.cy = parameters[1];
  distortion_.unproj_coeffs[0] = parameters[2];
  distortion_.unproj_coeffs[1] = parameters[3];
  distortion_.unproj_coeffs[2] = parameters[4];
  distortion_.unproj_coeffs[3] = parameters[5];
  distortion_.unproj_coeffs[4] = parameters[6];

  estimate_projection_coefficients();
}

void OcamLib::print() const
{
  std::cout << model_name_ << " parameters: "
            << "c=" << distortion_.c << ", "
            << "d=" << distortion_.d << ", "
            << "e=" << distortion_.e << ", "
            << "cx=" << common_params_.cx << ", "
            << "cy=" << common_params_.cy << std::endl;
  std::cout << "projection coefficients" << std::endl;
  for (auto elem : distortion_.proj_coeffs) {
    std::cout << elem << ", ";
  }
  std::cout << std::endl << "unprojection coefficients" << std::endl;
  for (auto elem : distortion_.unproj_coeffs) {
    std::cout << elem << ", ";
  }
  std::cout << std::endl;
}

void OcamLib::save_result(const std::string & result_path) const
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
  out << YAML::Key << "c" << YAML::Value << distortion_.c;
  out << YAML::Key << "d" << YAML::Value << distortion_.d;
  out << YAML::Key << "e" << YAML::Value << distortion_.e;
  out << YAML::Key << "cx" << YAML::Value << common_params_.cx;
  out << YAML::Key << "cy" << YAML::Value << common_params_.cy;

  out << YAML::Key << "coefficients" << YAML::Value;
  out << YAML::BeginMap;

  out << YAML::Key << "projection" << YAML::Value << YAML::Flow << YAML::BeginSeq;
  for (const auto & value : distortion_.proj_coeffs) {
    out << value;
  }
  out << YAML::EndSeq;

  out << YAML::Key << "unprojection" << YAML::Value << YAML::Flow << YAML::BeginSeq;
  for (const auto & value : distortion_.unproj_coeffs) {
    out << value;
  }
  out << YAML::EndSeq;

  out << YAML::EndMap;
  out << YAML::EndMap;

  out << YAML::EndMap;

  std::ofstream fout(result_path + "/" + model_name_ + ".yml");
  fout << out.c_str();
  fout << std::endl;
}

void OcamLib::estimate_projection_coefficients()
{
  constexpr auto MAX_POLY_NUM = 12;

  std::map<uint32_t, std::vector<double>> coefficient_candidates;
  double max_error = std::numeric_limits<double>::max();

  std::vector<Eigen::Vector3d> point3d_vec;
  std::vector<Eigen::Vector2d> point2d_vec;
  for (const auto & point2d : point2d_vec_) {
    const auto point3d = this->unproject(point2d, true);
    if (point3d.z() > 0.0) {
      point3d_vec.emplace_back(point3d);
      point2d_vec.emplace_back(point2d);
    }
  }

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
      const double theta = atan(-Z / r);

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
      b[2 * i] = u - common_params_.cx;
      b[(2 * i) + 1] = v - common_params_.cy;
    }
    const Eigen::VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
    distortion_.proj_coeffs.clear();
    distortion_.proj_coeffs = std::vector<double>(x.data(), x.data() + x.size());

    const double error = calculate_average_error(point3d_vec, point2d_vec);
    if (error < max_error) {
      max_error = error;
      coefficient_candidates.insert(std::make_pair(n, distortion_.proj_coeffs));
      best_poly_num = n;
    }
  }
  distortion_.proj_coeffs.clear();
  distortion_.proj_coeffs = coefficient_candidates.at(best_poly_num);
}

double OcamLib::calculate_average_error(
  const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
  double distance_sum = 0.0;
  for (auto i = 0U; i < point3d_vec.size(); ++i) {
    const auto & point3d = point3d_vec.at(i);
    const auto & point2d = point2d_vec.at(i);

    const Eigen::Vector2d estimated_point2d = project(point3d, true);
    distance_sum += (point2d - estimated_point2d).norm();
  }
  return distance_sum / point3d_vec.size();
}

void OcamLib::evaluate(const model::Base * const gt)
{
  const OcamLib * gt_model = dynamic_cast<const OcamLib *>(gt);

  const auto & est_pinhole_params = this->common_params_;
  const auto & gt_pinhole_params = gt_model->get_common_params();
  const auto & est_distortion_params = this->distortion_;
  const auto & gt_distortion_params = gt_model->get_distortion_params();

  const double diff_cx = est_pinhole_params.cx - gt_pinhole_params.cx;
  const double diff_cy = est_pinhole_params.cy - gt_pinhole_params.cy;

  const double diff_c = est_distortion_params.c - gt_distortion_params.c;
  const double diff_d = est_distortion_params.d - gt_distortion_params.d;
  const double diff_e = est_distortion_params.e - gt_distortion_params.e;

  double squared_norm_diff_distortion_coeffs = 0.0;
  for (auto i = 0; i < est_distortion_params.unproj_coeffs.size(); ++i) {
    const double diff =
      est_distortion_params.unproj_coeffs[i] - gt_distortion_params.unproj_coeffs[i];

    squared_norm_diff_distortion_coeffs += (diff * diff);
  }

  const double params_diff_norm = std::sqrt(
    diff_cx * diff_cx + diff_cy * diff_cy + diff_c * diff_c + diff_d * diff_d + diff_e * diff_e +
    squared_norm_diff_distortion_coeffs);

  std::cout << "parameter error: " << params_diff_norm << std::endl;
}
}  // namespace model
}  // namespace FCA
