cmake_minimum_required(VERSION 3.10)
set(PACKAGE_NAME fisheye-calib-adapter)
project(${PACKAGE_NAME})

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Ceres REQUIRED)
find_package(OpenCV REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

set(SOURCE_FILES
        src/adapter.cpp
        src/utils.cpp
        src/model/EUCM.cpp
        )
add_library(${PACKAGE_NAME}
        ${SOURCE_FILES}
        )

add_executable(main example/main.cpp)

target_link_libraries(main
    Eigen3::Eigen
    ${CERES_LIBRARIES}
    ${OpenCV_LIBS}
    ${PACKAGE_NAME}
)

target_include_directories(main PRIVATE
    ${CERES_INCLUDE_DIRS}
    ${OpenCV_INCLUDE_DIRS}
)
