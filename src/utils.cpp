#include "utils.hpp"

#include "model/EUCM.hpp"

namespace FCA
{
void Parse(
  int argc, char ** argv, std::string & input_model_name, std::string & output_model_name,
  std::string & config_path, std::string & dataset_path)
{
  // Default values
  input_model_name = "KB";
  output_model_name = "EUCM";
  dataset_path = "/config";
  dataset_path = "/dataset";

  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "-i" && i + 1 < argc) {
      input_model_name = argv[++i];
    } else if (std::string(argv[i]) == "-o" && i + 1 < argc) {
      output_model_name = argv[++i];
    } else if (std::string(argv[i]) == "-c" && i + 1 < argc) {
      config_path = argv[++i];
    } else if (std::string(argv[i]) == "-d" && i + 1 < argc) {
      dataset_path = argv[++i];
    }
  }
}

FisheyeCameraModelPtr Create(const std::string & model_name, const std::string & config_path)
{
  FisheyeCameraModelPtr model;
  if (model_name == "EUCM") {
    model = std::make_unique<model::EUCM>(model_name, config_path);
  } else if (model_name == "KB8") {
    // model = std::make_unique<model::KB8>(model_name, config_path);
  }
  return std::move(model);
}
}  // namespace FCA
