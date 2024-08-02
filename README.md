# Fisheye-Calib-Adapter

## Overview
The Fisheye-Calib-Adapter is a conversion tool designed to facilitate seamless transitions between different fisheye camera models. This tool is invaluable for applications in fields like robotics and autonomous driving, where fisheye cameras are increasingly common. It simplifies the process of camera model conversion, eliminating the need for cumbersome recalibrations. The adapter is user-friendly, straightforward, yet precise, offering conversion capabilities for a broader range of models compared to existing tools. It ensures that converted models function correctly in real-world applications such as SLAM, allowing researchers to easily translate input parameters into output parameters without the need for an image set or any recalibration processes.

## Paper

The pre-print paper related to this project can be found at: 

[Fisheye-Calib-Adapter: An Easy Tool for Fisheye Camera Model Conversion](https://arxiv.org/abs/2407.12405)

This work has been submitted to the IEEE for possible publication. Copyright may be transferred without notice, after which this version may no longer be accessible.

## Supported Camera Models
The Fisheye-Calib-Adapter supports a wide range of fisheye camera models, including:

- **UCM (Unified Camera Model)**
- **EUCM (Enhanced Unified Camera Model)**
- **DS (Double Sphere)**
- **KB (Kannala-Brandt)**
- **OCamCalib**
- **RadTan (Radial-Tangential Distortion Model)**

For those looking to implement a customized camera model, you can utilize the interface provided in '**include/model/custom.hpp**'. This allows for the extension and adaptation of the Fisheye-Calib-Adapter to meet specific requirements or to incorporate unique camera models not originally covered by the standard implementation.

## Requirements
The system requires the following libraries:

- Eigen
- OpenCV
- Ceres Solver
- yaml-cpp

## Setup
To set up the environment, follow these steps using Docker:

1. Build the Docker image:
   
```
docker build -t fca -f dockerfile/DockerFile .
```

3. Set display permissions for Docker and run the container:

```
xhost +local:docker
docker run -it -e DISPLAY=unix$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix -v $(pwd):/workspace fca
```

## Building the Project
To build the project, execute the following commands:
```
mkdir build
cd build
cmake ..
make -j 4
```
## Running the Tool
To run the tool, use the following command:
```
./main example/config.yml
```

You can figure out the input-output model in the example/config.yml described as below

## Configuration File Description
The config.yml file should include the following parameters:

- **dataset_path**: The root path where the fisheye model is stored.
- **input_model**, output_model: The input and output models you wish to convert.
- **image_name**: The name of the image to apply if you want, which must be located within dataset_path.
- **show_image**: Set to true if you want to display the image.
- **sample_point**: The number of point samples.
- **result_path**: The path where the output model will be stored.
- **save_result**: Set to true to save the output model.
- **display_optimization_progress**: Displays the optimization process if set to true.
- **evaluation**: Set to true to perform evaluation (calculates reprojection error, parameter error, RMSE/SSIM if an image is set).

## Citation

If you use Fisheye-Calib-Adapter, please use the following BibTeX entry.

```
@article{FCA2024,
  title={Fisheye-Calib-Adapter: An Easy Tool for Fisheye Camera Model Conversion},
  author={Sangjun Lee},
  journal={arXiv preprint},
  year={2024}
}
```
