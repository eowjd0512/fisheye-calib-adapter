// #include <ceres/ceres.h>

// #include <Eigen/Dense>
// #include <iostream>
// #include <opencv2/opencv.hpp>

#include <string>

#include "adapter.hpp"
#include "model/base.hpp"
#include "utils.hpp"

int main(int argc, char ** argv)
{
  std::string input_model_name;
  std::string output_model_name;
  std::string config_path;
  std::string dataset_path;

  FCA::Parse(argc, argv, input_model_name, output_model_name, config_path, dataset_path);

  FCA::FisheyeCameraModelPtr input_model = FCA::Create(input_model_name, config_path);
  FCA::FisheyeCameraModelPtr output_model = FCA::Create(output_model_name, config_path);

  FCA::Adapter adapter(dataset_path);

  bool show = true;
  adapter.adapt(*input_model, *output_model, show);

  output_model->print();
  return 0;
}
