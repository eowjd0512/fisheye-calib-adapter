#ifndef MODEL__BASE_HPP_
#define MODEL__BASE_HPP_

#include <yaml-cpp/yaml.h>

#include <eigen3/Eigen/Dense>
#include <iostream>
#include <memory>
#include <string>

namespace FCA
{
namespace model
{
class Base
{
public:
  Base(const std::string & model_name, const std::string & config_path)
  : model_name_(model_name), config_path_(config_path)
  {
  }

  virtual ~Base() {}

  virtual void parse() = 0;
  virtual void initialize() = 0;
  virtual Eigen::Vector2d project(const Eigen::Vector3d & point3d) const = 0;
  virtual Eigen::Vector3d unproject(const Eigen::Vector2d & point2d) const = 0;
  virtual Eigen::MatrixXd calculate_jacobian() const = 0;
  virtual void optimize() = 0;
  virtual void print() const = 0;

protected:
  std::string model_name_;
  std::string config_path_;
};
}  // namespace model
using FisheyeCameraModel = model::Base;
using FisheyeCameraModelPtr = std::unique_ptr<model::Base>;
}  // namespace FCA
#endif
