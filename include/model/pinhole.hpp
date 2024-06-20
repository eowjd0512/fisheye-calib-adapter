#ifndef MODEL__PINHOLE_HPP_
#define MODEL__PINHOLE_HPP_

#include <string>

#include "model/base.hpp"

namespace FCA
{
namespace model
{
class Pinhole : public Base
{
public:

  Pinhole(const std::string & model_name, const std::string & config_path);

  void parse() override;
  void set_sample_points(const std::vector<Eigen::Vector2d> & point2d_vec) override;
  void initialize(
    const Base::Params & common_params, const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec) override;
  Eigen::Vector2d project(const Eigen::Vector3d & point3d) const override;
  Eigen::Vector3d unproject(const Eigen::Vector2d & point2d) const override;
  void optimize(
    const std::vector<Eigen::Vector3d> & point3d_vec,
    const std::vector<Eigen::Vector2d> & point2d_vec) override;
  void print() const override;
  void save_result(const std::string& result_path) const override;

private:
};

}  // namespace model
}  // namespace FCA
#endif
