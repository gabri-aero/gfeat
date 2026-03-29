cmake_minimum_required(VERSION 3.21)
project(gfeatpy VERSION ${PROJECT_VERSION} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Root directories
set(PROJECT_ROOT ${CMAKE_SOURCE_DIR})

# Data directory
set(DATA_DIR "${PROJECT_ROOT}/data")
add_definitions(-DDATA_DIR="${DATA_DIR}")

# Check if building Python wrapper or C++ tests
if(BUILD_GFEATPY)
    message(STATUS "Building python library")
    set(BUILD_CPP_TESTS OFF)
else()
    message(STATUS "Building C++ tests")
    set(BUILD_CPP_TESTS ON)
endif()

# ---------------------------
# Header-only gfeat library
# ---------------------------
set(CPP_TARGET gfeat)
file(GLOB_RECURSE HEADER_FILES "${PROJECT_ROOT}/include/gfeat/*.hpp")
add_library(${CPP_TARGET} INTERFACE gfeat ${HEADER_FILES})

# Compile options
if(MSVC)
    target_compile_definitions(${CPP_TARGET} INTERFACE _USE_MATH_DEFINES)
    target_compile_options(${CPP_TARGET} INTERFACE /permissive- /O2 /fp:fast)
else()
    target_compile_options(${CPP_TARGET} INTERFACE -O3 -ffast-math)
endif()

# ---------------------------
# C++ tests
# ---------------------------
if(BUILD_CPP_TESTS)
    enable_testing()
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(googletest)
    include(GoogleTest)
    add_subdirectory(tests)
endif()

# ---------------------------
# Python wrapper _core
# ---------------------------
if(BUILD_GFEATPY)
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
    set(PYBIND11_FINDPYTHON ON)
    find_package(pybind11 CONFIG REQUIRED)
    set(PYTHON_SITE_PACKAGES_DIR ${Python3_SITE_PACKAGES})

    # Define the Python library generated from C++
    python_add_library(_core MODULE gfeatpy/main.cpp WITH_SOABI)

    # Add dependency on the header-only library
    add_dependencies(_core ${CPP_TARGET})

    # Link the header-only library to the Python module
    target_link_libraries(_core PRIVATE ${CPP_TARGET})

    # Settings
    target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})
    set_target_properties(_core PROPERTIES
    INSTALL_RPATH "$ORIGIN"
    BUILD_WITH_INSTALL_RPATH TRUE
    )

    # Link root directory
    target_include_directories(_core PRIVATE ${PROJECT_ROOT})

    # Link pybind11 headers
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
    file(GLOB PYTHON_FILES "${PROJECT_ROOT}/gfeatpy/*.py")
    install(FILES ${PYTHON_FILES} DESTINATION gfeatpy)

endif()