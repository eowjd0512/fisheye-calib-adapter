/*
MIT License

Copyright (c) 2024 Sangjun Lee, STRADVISION

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <chrono>
#include <filesystem>
#include <opencv2/core/utils/logger.hpp>
#include <string>

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
  FCA::FisheyeCameraModelPtr gt_output_model =
    FCA::Create(output_model_name, dataset_path.string());

  FCA::Adapter::Config fca_config;
  fca_config.sample_point = config["sample_point"].as<int32_t>();
  // fca_config.sample_point = i;
  fca_config.display_optimization_progress = config["display_optimization_progress"].as<bool>();
  FCA::Adapter adapter(fca_config, input_model.get(), output_model.get());

  const std::string image_name = config["image_name"].as<std::string>();
  adapter.set_image((dataset_path / image_name).string());

  // Adapt!
  std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
  adapter.adapt();
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  auto sec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "computing time(ms): " << sec.count() << std::endl;

  if (config["show_image"].as<bool>()) {
    adapter.show_image(result_path.string());
  }

  if (config["save_result"].as<bool>()) {
    output_model->save_result(result_path.string());
  }

  if (config["evaluation"].as<bool>()) {
    adapter.evaluate(gt_output_model.get());
  }

  return 0;
}
