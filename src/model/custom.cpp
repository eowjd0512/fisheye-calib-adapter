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

#include "model/custom.hpp"

namespace FCA
{
namespace model
{

Custom::Custom(const std::string & model_name, const std::string & config_path)
: Base(model_name, config_path)
{
}

void Custom::parse() {}
void Custom::set_sample_points(const std::vector<Eigen::Vector2d> & point2d_vec){};
void Custom::initialize(
  const Base::Params & common_params, const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec)
{
  assert(point3d_vec.size() == point2d_vec.size());
}
Eigen::Vector2d Custom::project(const Eigen::Vector3d & point3d, bool condition) const { return Eigen::Vector2d(); }
Eigen::Vector3d Custom::unproject(const Eigen::Vector2d & point2d, bool condition) const
{
  return Eigen::Vector3d();
}
void Custom::optimize(
  const std::vector<Eigen::Vector3d> & point3d_vec,
  const std::vector<Eigen::Vector2d> & point2d_vec,
  bool display_optimization_progress)
{
}
void Custom::print() const {}

void Custom::save_result(const std::string& result_path) const {}

void Custom::evaluate(const model::Base * const gt){}
}  // namespace model
}  // namespace FCA
