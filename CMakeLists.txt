cmake_minimum_required(VERSION 3.21)
project(gfeatpy LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the root source directory
set(PROJECT_ROOT ${CMAKE_SOURCE_DIR})

# Define data directory path
set(DATA_DIR "${PROJECT_ROOT}/data")
add_definitions(-DDATA_DIR="${DATA_DIR}")

# Check if C++ tests should be built
if(BUILD_GFEATPY)
message(STATUS "Building python library")
set(BUILD_CPP_TESTS OFF)
else()
message(STATUS "Building C++ tests")
set(BUILD_CPP_TESTS ON)
endif()

# Retrieve all source files
file(GLOB_RECURSE SOURCES "include/*.hpp")
message(STATUS "Files found: ${SOURCES}")
# Create C++ library
set(CPP_TARGET gfeat)
add_library(${CPP_TARGET} SHARED ${SOURCES})
set_target_properties(${CPP_TARGET} PROPERTIES LINKER_LANGUAGE CXX)
add_compile_options(-march=native -ffast-math)  # -O3 might be unsuitable for debugging
target_include_directories(${CPP_TARGET} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen)

# Setup external libraries

if(BUILD_CPP_TESTS)
enable_testing()
# Google Test
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(googletest)
target_link_libraries(${CPP_TARGET} PRIVATE GTest::gtest_main)
include(GoogleTest)
# Create C++ tests
add_subdirectory(tests)
endif()

# Build python wrapper _core
if(BUILD_GFEATPY)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)
set(PYTHON_SITE_PACKAGES_DIR ${Python3_SITE_PACKAGES})

# Define the Python library generated from C++
python_add_library(_core MODULE gfeatpy/main.cpp WITH_SOABI)

# Settings
target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})
target_compile_options(_core PRIVATE -O3 -march=native -ffast-math -fopenmp)
set_target_properties(_core PROPERTIES
  INSTALL_RPATH "$ORIGIN"
  BUILD_WITH_INSTALL_RPATH TRUE
)

# Link 
target_include_directories(_core PRIVATE ${PROJECT_ROOT})
target_link_libraries(_core PRIVATE pybind11::headers)

# OpenMP setup
find_package(OpenMP REQUIRED)
target_link_libraries(_core PUBLIC OpenMP::OpenMP_CXX)

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

endif()