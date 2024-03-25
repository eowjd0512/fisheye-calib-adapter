#ifndef ADAPTER_HPP_
#define ADAPTER_HPP_

#include <string>

#include "model/base.hpp"

namespace FCA
{
class Adapter
{
public:
  Adapter(
    const FCA::FisheyeCameraModel * const input_model,
    FCA::FisheyeCameraModel * const output_model);

  void adapt();

  void compare_image(const std::string & image_path, bool show = false);

private:
  const FCA::FisheyeCameraModel * const input_model_;
  FCA::FisheyeCameraModel * const output_model_;
};
}  // namespace FCA
#endif
