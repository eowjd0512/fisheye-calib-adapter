#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <string>

#include "model/base.hpp"

namespace FCA
{
void Parse(
  int argc, char ** argv, std::string & input_model_name, std::string & output_model_name,
  std::string & dataset_path, std::string & result_path);

FisheyeCameraModelPtr Create(const std::string & model_name, const std::string & dataset_path);

}  // namespace FCA
#endif
