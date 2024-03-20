#include "adapter.hpp"

namespace FCA
{

Adapter::Adapter(const std::string & dataset_path) : dataset_path_(dataset_path) {}

void Adapter::adapt(
  const FisheyeCameraModel & input_model, const FisheyeCameraModel & output_model, bool show)
{
}

}  // namespace FCA