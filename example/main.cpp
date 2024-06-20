#include <filesystem>
#include <string>

#include <opencv2/core/utils/logger.hpp>

#include "adapter.hpp"
#include "model/base.hpp"
#include "utils.hpp"

namespace fs = std::filesystem;

int main(int argc, char ** argv)
{
   cv::utils::logging::setLogLevel(cv::utils::logging::LOG_LEVEL_SILENT);

  fs::path project_dir = PROJECT_SOURCE_DIR;
  fs::path config_path = project_dir / argv[1];

  const YAML::Node config = YAML::LoadFile(config_path.string());

  fs::path dataset_path = project_dir / config["dataset_path"].as<std::string>();
  fs::path result_path = project_dir / config["result_path"].as<std::string>();

  const std::string input_model_name = config["input_model"].as<std::string>();
  const std::string output_model_name = config["output_model"].as<std::string>();

  FCA::FisheyeCameraModelPtr input_model = FCA::Create(input_model_name, dataset_path.string());
  FCA::FisheyeCameraModelPtr output_model = FCA::Create(output_model_name, dataset_path.string());

  FCA::Adapter::Config fca_config;
  fca_config.sample_point = config["sample_point"].as<int32_t>();
  fca_config.display_optimization_progress = config["display_optimization_progress"].as<bool>();
  FCA::Adapter adapter(fca_config, input_model.get(), output_model.get());

  const std::string image_name = config["image_name"].as<std::string>();
  adapter.set_image((dataset_path / image_name).string());

  // Adapt!
  adapter.adapt();

  if (config["show_image"].as<bool>()) {
    adapter.show_image();
  }

  if (config["save_result"].as<bool>()) {
    output_model->save_result(result_path.string());
  }

  if (config["evaluation"].as<bool>()) {
    adapter.evaluate();
  }
  return 0;
}
