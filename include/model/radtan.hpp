#ifndef MODEL__RADTAN_HPP_
#define MODEL__RADTAN_HPP_

#include <string>

#include "model/base.hpp"

namespace FCA
{
namespace model
{
class Radtan : public Base
{
public:
  struct Params
  {
    double k1;
    double k2;
    double k3;
    double p1;
    double p2;
  };

  Radtan(const std::string & model_name, const std::string & config_path);

  void parse() override;
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
  Params distortion_;
};

}  // namespace model
}  // namespace FCA
#endif
