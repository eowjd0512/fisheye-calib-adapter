#ifndef MODEL__BASE_HPP_
#define MODEL__BASE_HPP_

#include <ceres/ceres.h>
#include <yaml-cpp/yaml.h>

#include <Eigen/Dense>
#include <fstream>
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
  struct Params
  {
    int32_t width;
    int32_t height;
    double fx;
    double fy;
    double cx;
    double cy;
  };

  Base(const std::string & model_name, const std::string & config_path)
  : model_name_(model_name), config_path_(config_path)
  {
  }

  virtual ~Base() {}

  virtual void parse() = 0;
  virtual void initialize(
    const Params & common_params, const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec) = 0;
  virtual Eigen::Vector2d project(const Eigen::Vector3d & point3d) const = 0;
  virtual Eigen::Vector3d unproject(const Eigen::Vector2d & point2d) const = 0;
  virtual void optimize(
    const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec) = 0;
  virtual void print() const = 0;
  virtual void save_result(const std::string & result_path) const = 0;

  const Params & get_common_params() const { return common_params_; };

protected:
  std::string model_name_;
  std::string config_path_;
  Params common_params_;
};
}  // namespace model
using FisheyeCameraModel = model::Base;
using FisheyeCameraModelPtr = std::unique_ptr<model::Base>;
}  // namespace FCA
#endif
