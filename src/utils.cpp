#include "utils.hpp"

#include "model/EUCM.hpp"
#include "model/KB8.hpp"
#include "model/OcamLib.hpp"
#include "model/UCM.hpp"
#include "model/double_sphere.hpp"

namespace FCA
{
void Parse(
  int argc, char ** argv, std::string & input_model_name, std::string & output_model_name,
  std::string & dataset_path, std::string & result_path)
{
  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "-i" && i + 1 < argc) {
      input_model_name = argv[++i];
    } else if (std::string(argv[i]) == "-o" && i + 1 < argc) {
      output_model_name = argv[++i];
    } else if (std::string(argv[i]) == "-d" && i + 1 < argc) {
      dataset_path = argv[++i];
    } else if (std::string(argv[i]) == "-r" && i + 1 < argc) {
      result_path = argv[++i];
    }
  }
}

FisheyeCameraModelPtr Create(const std::string & model_name, const std::string & dataset_path)
{
  FisheyeCameraModelPtr model;

  if (model_name == "EUCM") {
    model = std::make_unique<model::EUCM>(model_name, dataset_path);
  } else if (model_name == "UCM") {
    model = std::make_unique<model::UCM>(model_name, dataset_path);
  } else if (model_name == "KB8") {
    model = std::make_unique<model::KB8>(model_name, dataset_path);
  } else if (model_name == "OcamLib") {
    model = std::make_unique<model::OcamLib>(model_name, dataset_path);
  } else if (model_name == "DS") {
    model = std::make_unique<model::DoubleSphere>(model_name, dataset_path);
  } else {
    std::cerr << "Use Fisheye Camera Model within {UCM, EUCM, DS, KB8, OcamLib}" << std::endl;
  }

  return std::move(model);
}
}  // namespace FCA
