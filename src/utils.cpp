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
