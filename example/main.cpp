#include <string>

#include "adapter.hpp"
#include "model/base.hpp"
#include "utils.hpp"

int main(int argc, char ** argv)
{
  // Default values
  std::string input_model_name = "KB8";
  std::string output_model_name = "DS";
  std::string dataset_path = std::string(PROJECT_SOURCE_DIR) + "/dataset";
  std::string result_path = std::string(PROJECT_SOURCE_DIR) + "/result";

  FCA::Parse(argc, argv, input_model_name, output_model_name, dataset_path, result_path);

  FCA::FisheyeCameraModelPtr input_model = FCA::Create(input_model_name, dataset_path);
  FCA::FisheyeCameraModelPtr output_model = FCA::Create(output_model_name, dataset_path);

  FCA::Adapter adapter(input_model.get(), output_model.get());
  adapter.adapt();
  adapter.evaluate();
  bool show = true;
  // adapter.compare_image(dataset_path + "/image.jpg", show);

  // output_model->save_result(result_path);

  return 0;
}

// #include <opencv2/opencv.hpp>
// #include <string>
// #include <iomanip>
// #include <sstream>

// int main() {
//     // 비디오 파일 경로
//     std::string dataset_path = std::string(PROJECT_SOURCE_DIR) + "/dataset/";
//     std::string videoPath = dataset_path + "sample_video.mp4";

//     // 비디오를 열기 위한 객체 생성
//     cv::VideoCapture cap(videoPath);

//     // 비디오 파일 열기 실패 시
//     if (!cap.isOpened()) {
//         std::cerr << "Error opening video file" << std::endl;
//         return -1;
//     }

//     cv::Mat frame;
//     int frameNumber = 0;

//     while (true) {
//         // 비디오에서 프레임을 읽음
//         cap >> frame;

//         // 프레임이 비어있으면 종료 (비디오 끝)
//         if (frame.empty()) {
//             break;
//         }

//         // 파일 이름을 000.jpg, 001.jpg 형식으로 생성
//         std::stringstream ss;
//         ss << std::setw(3) << std::setfill('0') << frameNumber << ".jpg";

//         // 이미지 파일로 저장
//         cv::imwrite(dataset_path + "output/" + ss.str(), frame);

//         // 프레임 번호 증가
//         frameNumber++;
//     }

//     return 0;
// }