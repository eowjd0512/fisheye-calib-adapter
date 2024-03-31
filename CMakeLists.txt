cmake_minimum_required(VERSION 3.10)
set(PACKAGE_NAME fisheye-calib-adapter)
project(${PACKAGE_NAME})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 3.3 REQUIRED)
find_package(Ceres REQUIRED)
find_package(OpenCV REQUIRED)
find_package(yaml-cpp REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${EIGEN3_INCLUDE_DIR} ${OpenCV_INCLUDE_DIRS})

set(SOURCE_FILES
        src/adapter.cpp
        src/utils.cpp
        src/model/double_sphere.cpp
        src/model/EUCM.cpp
        src/model/KB8.cpp
        src/model/OcamLib.cpp
        src/model/UCM.cpp
        )
add_library(${PACKAGE_NAME}
        ${SOURCE_FILES}
        )

add_definitions(-DPROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

add_executable(main example/main.cpp)

target_link_libraries(main
    ${PACKAGE_NAME}
    Eigen3::Eigen
    ${CERES_LIBRARIES}
    ${OpenCV_LIBS}
    ${YAML_CPP_LIBRARIES}
)

target_include_directories(main PRIVATE
    ${CERES_INCLUDE_DIRS}
    ${OpenCV_INCLUDE_DIRS}
    ${YAML_CPP_INCLUDE_DIR}
)
