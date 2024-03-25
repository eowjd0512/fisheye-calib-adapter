
// #include <opencv2/opencv.hpp>

#include <string>

#include "adapter.hpp"
#include "model/base.hpp"
#include "utils.hpp"

int main(int argc, char ** argv)
{
  // Default values
  std::string input_model_name = "OcamLib";  //"KB8";
  std::string output_model_name = "EUCM";
  std::string dataset_path = std::string(PROJECT_SOURCE_DIR) + "/dataset";

  FCA::Parse(argc, argv, input_model_name, output_model_name, dataset_path);

  FCA::FisheyeCameraModelPtr input_model = FCA::Create(input_model_name, dataset_path);
  FCA::FisheyeCameraModelPtr output_model = FCA::Create(output_model_name, dataset_path);

  FCA::Adapter adapter(input_model.get(), output_model.get());
  adapter.adapt();

  bool show =  true;
  adapter.compare_image(dataset_path + "/image.jpg", show);

  output_model->print();
  return 0;
}
