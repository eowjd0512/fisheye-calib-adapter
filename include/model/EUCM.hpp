#ifndef MODEL__EUCM_HPP_
#define MODEL__EUCM_HPP_

#include <string>

#include "model/base.hpp"

namespace FCA
{
namespace model
{
class EUCM : public Base
{
public:
  EUCM(const std::string & model_name, const std::string & config_path);

  void parse() override;
  void initialize() override;
  void optimize() override;
  Eigen::Vector2d project() const override;
  Eigen::Vector3d unproject() const override;
  Eigen::MatrixXd calculate_jacobian() const override;
  void print() const override;
};
}  // namespace model
}  // namespace FCA
#endif
