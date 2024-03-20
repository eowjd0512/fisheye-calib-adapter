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

  void print() const override;
  void project() const override;
  void unproject() const override;
  void calculate_jacobian() const override;
  void get_jacobian() const override;
};
}  // namespace model
using FisheyeCameraModel = model::Base;
}  // namespace FCA
#endif
