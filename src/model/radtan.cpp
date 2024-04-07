#include "model/radtan.hpp"

namespace FCA
{
namespace model
{

Radtan::Radtan(const std::string & model_name, const std::string & config_path)
: Base(model_name, config_path)
{
  parse();
}

void Radtan::parse()
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
  distortion_.k1 = config["parameter"]["k1"].as<double>();
  distortion_.k2 = config["parameter"]["k2"].as<double>();
  distortion_.k3 = config["parameter"]["k3"].as<double>();
  distortion_.p1 = config["parameter"]["p1"].as<double>();
  distortion_.p2 = config["parameter"]["p2"].as<double>();
}

void Radtan::initialize(
  const Base::Params & common_params, const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
  assert(point3d_vec.size() == point2d_vec.size());
  std::cerr << ERROR_STR("NOT SUPPORTED YET") << std::endl;
}

Eigen::Vector2d Radtan::project(const Eigen::Vector3d & point3d) const
{
  if (point3d.z() <= 0.0) {
    return {-1., -1.};
  }

  const double x = point3d.x() / point3d.z();
  const double y = point3d.y() / point3d.z();

  const double r2 = x * x + y * y;
  const double radial_distortion =
    1.0 + r2 * (distortion_.k1 + r2 * (distortion_.k2 + r2 * distortion_.k3));
  const double x_tangential = 2.0 * distortion_.p1 * x * y + distortion_.p2 * (r2 + 2.0 * x * x);
  const double y_tangential = distortion_.p1 * (r2 + 2.0 * y * y) + 2.0 * distortion_.p2 * x * y;
  const double x_distorted = x * radial_distortion + x_tangential;
  const double y_distorted = y * radial_distortion + y_tangential;

  const double u = common_params_.fx * x_distorted + common_params_.cx;
  const double v = common_params_.fy * y_distorted + common_params_.cy;

  return {u, v};
}

Eigen::Vector3d Radtan::unproject(const Eigen::Vector2d & point2d) const
{
  const double x = (point2d.x() - common_params_.cx) / common_params_.fx;
  const double y = (point2d.y() - common_params_.cy) / common_params_.fy;

  const double r2 = x * x + y * y;
  const double radial_distortion =
    1.0 + r2 * (distortion_.k1 + r2 * (distortion_.k2 + r2 * distortion_.k3));
  const double x_tangential = 2.0 * distortion_.p1 * x * y + distortion_.p2 * (r2 + 2.0 * x * x);
  const double y_tangential = distortion_.p1 * (r2 + 2.0 * y * y) + 2.0 * distortion_.p2 * x * y;
  const double x_corrected = x * radial_distortion + x_tangential;
  const double y_corrected = y * radial_distortion + y_tangential;

  const double inv_norm = 1.0 / sqrt(x_corrected * x_corrected + y_corrected * y_corrected + 1.0);

  return {x_corrected * inv_norm, y_corrected * inv_norm, inv_norm};
}

void Radtan::optimize(
  const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
  std::cerr << ERROR_STR("NOT SUPPORTED YET") << std::endl;
}

void Radtan::print() const
{
  std::cout << "Final parameters: "
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

void Radtan::save_result(const std::string & result_path) const
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
}  // namespace model
}  // namespace FCA
