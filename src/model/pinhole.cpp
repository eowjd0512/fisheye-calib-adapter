#include "model/pinhole.hpp"

namespace FCA
{
namespace model
{
Pinhole::Pinhole(const std::string & model_name, const std::string & config_path)
: Base(model_name, config_path)
{
}

void Pinhole::parse()
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
}

void Pinhole::set_sample_points(const std::vector<Eigen::Vector2d> & point2d_vec){};

void Pinhole::initialize(
  const Base::Params & common_params, const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
  assert(point3d_vec.size() == point2d_vec.size());
  std::cerr << ERROR_STR("NOT SUPPORTED YET") << std::endl;
}

Eigen::Vector2d Pinhole::project(const Eigen::Vector3d & point3d, bool condition) const
{
  if (point3d.z() <= 0.0) {
    return {-1., -1.};
  }

  const double u = common_params_.fx * point3d.x() / point3d.z() + common_params_.cx;
  const double v = common_params_.fy * point3d.y() / point3d.z() + common_params_.cy;
  return {u, v};
}

Eigen::Vector3d Pinhole::unproject(const Eigen::Vector2d & point2d, bool condition) const
{
  const double x = (point2d.x() - common_params_.cx) / common_params_.fx;
  const double y = (point2d.y() - common_params_.cy) / common_params_.fy;

  Eigen::Vector3d ray(x, y, 1.0);
  ray.normalize();

  return ray;
}

void Pinhole::optimize(
  const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec, bool display_optimization_progress)
{
  std::cerr << ERROR_STR("NOT SUPPORTED YET") << std::endl;
}

void Pinhole::print() const
{
  std::cout << model_name_ << " parameters: "
            << "fx=" << common_params_.fx << ", "
            << "fy=" << common_params_.fy << ", "
            << "cx=" << common_params_.cx << ", "
            << "cy=" << common_params_.cy << std::endl;
}

void Pinhole::save_result(const std::string & result_path) const
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
  out << YAML::EndMap;

  out << YAML::EndMap;

  std::ofstream fout(result_path + "/" + model_name_ + ".yml");
  fout << out.c_str();
  fout << std::endl;
}

void Pinhole::evaluate(const model::Base * const gt) {}
}  // namespace model
}  // namespace FCA
