#ifndef MODEL__BASE_HPP_
#define MODEL__BASE_HPP_

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

  virtual void print() const = 0;
  virtual void project() const = 0;
  virtual void unproject() const = 0;
  virtual void calculate_jacobian() const = 0;
  virtual void get_jacobian() const = 0;

protected:
  std::string model_name_;
  std::string config_path_;
};
}  // namespace model
using FisheyeCameraModel = model::Base;
using FisheyeCameraModelPtr = std::unique_ptr<model::Base>;
}  // namespace FCA
#endif
