#ifndef ADAPTER_HPP_
#define ADAPTER_HPP_

#include <string>

#include "model/base.hpp"

namespace FCA
{
class Adapter
{
public:
  Adapter(const std::string & dataset_path);

  void adapt(
    const FisheyeCameraModel & input_model, const FisheyeCameraModel & output_model,
    bool show = false);

private:
  std::string dataset_path_;
};
}  // namespace FCA
#endif
