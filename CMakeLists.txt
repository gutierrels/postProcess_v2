
cmake_minimum_required(VERSION 3.16)
project(postProcess LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(penred2Coincidence
    src/main.cpp
    src/config/config.cpp
    src/geometry/geometry.cpp
    src/geometry/projection.cpp
    src/io/singles_reader.cpp
    src/io/pair_list.cpp
    src/coincidence/coincidence_engine.cpp
    src/histogram/histograms.cpp
    src/histogram/histogram_writer.cpp
    src/output/lm_writer.cpp
)

target_include_directories(penred2Coincidence PRIVATE src)

# Optional: Add release optimization flags if compiling in Release mode
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(penred2Coincidence PRIVATE -O3 -march=native)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(penred2Coincidence PRIVATE -g -O0 -fno-omit-frame-pointer -fsanitize=address,undefined)
    target_link_options(penred2Coincidence PRIVATE -fsanitize=address,undefined)
endif()
