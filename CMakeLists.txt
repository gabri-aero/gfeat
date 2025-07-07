cmake_minimum_required(VERSION 3.21)
project(gfeatpy LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the root source directory
set(PROJECT_ROOT ${CMAKE_SOURCE_DIR})

# Define data directory path
set(DATA_DIR "${PROJECT_ROOT}/data")
add_definitions(-DDATA_DIR="${DATA_DIR}")

# Find required packages
find_package(pybind11 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Python3 REQUIRED)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
set(PYTHON_SITE_PACKAGES_DIR ${Python3_SITE_PACKAGES})

# Define the Python library generated from C++
python_add_library(_core MODULE gfeatpy/main.cpp WITH_SOABI)

# Settings
target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})
set_target_properties(_core PROPERTIES
  INSTALL_RPATH "$ORIGIN"
  BUILD_WITH_INSTALL_RPATH TRUE
)

# Link 
target_include_directories(_core PRIVATE ${PROJECT_ROOT})
target_include_directories(_core PRIVATE ${EIGEN3_INCLUDE_DIR})
target_link_libraries(_core PRIVATE pybind11::headers)

# Link external directories
file(GLOB EXTERNAL_DIRS LIST_DIRECTORIES true "${PROJECT_ROOT}/external/*")
foreach(DIR ${EXTERNAL_DIRS})
    if(IS_DIRECTORY ${DIR})
        message(STATUS "Adding external directory: ${DIR}")
        target_include_directories(_core PRIVATE ${DIR})
    endif()
endforeach()

# Install targets
install(TARGETS _core DESTINATION gfeatpy)
file(GLOB PYTHON_FILES "${PROJECT_ROOT}/py/*.py")
install(FILES ${PYTHON_FILES} DESTINATION gfeatpy)