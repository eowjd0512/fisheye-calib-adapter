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
  YAML::Node config = YAML::LoadFile(config_path_ + "/" + model_name_ + ".yml");

  // Read image dimensions
  int32_t width = config["image"]["width"].as<int32_t>();
  int32_t height = config["image"]["height"].as<int32_t>();

  // Read coefficients
  std::vector<double> projection_coefficients =
    config["coefficients"]["projection"].as<std::vector<double>>();
  std::vector<double> unprojection_coefficients =
    config["coefficients"]["unprojection"].as<std::vector<double>>();

  // Read parameters
  std::vector<double> affine_parameters = config["parameter"]["affine"].as<std::vector<double>>();
  int32_t cx = config["parameter"]["cx"].as<int32_t>();
  int32_t cy = config["parameter"]["cy"].as<int32_t>();
}
void OcamLib::initialize() {}
void OcamLib::optimize() {}
Eigen::Vector2d OcamLib::project() const {}
Eigen::Vector3d OcamLib::unproject() const {}
Eigen::MatrixXd OcamLib::calculate_jacobian() const {}
void OcamLib::print() const {}
}  // namespace model
}  // namespace FCA
