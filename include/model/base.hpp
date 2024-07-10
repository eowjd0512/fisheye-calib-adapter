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

#ifndef MODEL__BASE_HPP_
#define MODEL__BASE_HPP_

#include <ceres/ceres.h>
#include <yaml-cpp/yaml.h>

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#define ERROR_STR(MSG)                                                                     \
  std::string(__FILE__) + ": " + std::to_string(__LINE__) + ", " + std::string(__func__) + \
    "(): " + MSG

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
  virtual void set_sample_points(const std::vector<Eigen::Vector2d> & point2d_vec) = 0;
  virtual void initialize(
    const Params & common_params, const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec) = 0;
  virtual Eigen::Vector2d project(const Eigen::Vector3d & point3d, bool condition = true) const = 0;
  virtual Eigen::Vector3d unproject(
    const Eigen::Vector2d & point2d, bool condition = true) const = 0;
  virtual void optimize(
    const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec,
    bool display_optimization_progress = false) = 0;
  virtual void print() const = 0;
  virtual void save_result(const std::string & result_path) const = 0;
  virtual void evaluate(const model::Base * const gt) = 0;
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
