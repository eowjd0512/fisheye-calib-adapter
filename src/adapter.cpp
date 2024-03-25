#include "adapter.hpp"

namespace FCA
{

Adapter::Adapter(
  const FCA::FisheyeCameraModel * const input_model, FCA::FisheyeCameraModel * const output_model)
: input_model_(input_model), output_model_(output_model)
{
}

void Adapter::adapt()
{
  // Sample points

  // Unproject from input model

  // associate <3D, 2D>

  // Initialize output model

  // Optimize output model
}

void Adapter::compare_image(const std::string & image_path, bool show)
{
  // if(!output_model_->is_estimated()){
  //   this->adapt();
  // }
}
}  // namespace FCA