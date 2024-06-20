#include "utils.hpp"

#include "model/EUCM.hpp"
#include "model/KB.hpp"
#include "model/OcamLib.hpp"
#include "model/UCM.hpp"
#include "model/double_sphere.hpp"

namespace FCA
{

FisheyeCameraModelPtr Create(const std::string & model_name, const std::string & dataset_path)
{
  FisheyeCameraModelPtr model;

  if (model_name == "EUCM") {
    model = std::make_unique<model::EUCM>(model_name, dataset_path);
  } else if (model_name == "UCM") {
    model = std::make_unique<model::UCM>(model_name, dataset_path);
  } else if (model_name == "KB") {
    model = std::make_unique<model::KB>(model_name, dataset_path);
  } else if (model_name == "OcamLib") {
    model = std::make_unique<model::OcamLib>(model_name, dataset_path);
  } else if (model_name == "DS") {
    model = std::make_unique<model::DoubleSphere>(model_name, dataset_path);
  } else {
    std::cerr << "Use Fisheye Camera Model within {UCM, EUCM, DS, KB, OcamLib}" << std::endl;
  }

  return std::move(model);
}
}  // namespace FCA
