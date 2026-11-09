cmake_minimum_required(VERSION 3.16...4.4)
project(jpegultrascancpp VERSION 0.1.0 LANGUAGES C CXX)

if(MSVC)
    # force static msvc runtime linkage (/mt for release, /mtd for debug)
    if(POLICY CMP0091)
        cmake_policy(SET CMP0091 NEW)
        set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif()
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# without this a plain `cmake -S . -B build` compiles with no -O flag at all
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
    message(STATUS "CMAKE_BUILD_TYPE not set, defaulting to Release")
endif()

option(JPEGULTRASCAN_BUILD_CLI "Build the command-line tool" ON)

include(cmake/libjpeg.cmake)

# std::thread needs an explicit -pthread on some platforms; freebsd fails to link without it
find_package(Threads REQUIRED)

add_library(jpegultrascan_core STATIC
    src/engine_direct.cpp
    src/jpeg_io.cpp
    src/jpegultrascan.cpp
    src/scan_script.cpp
    src/search.cpp
)
target_include_directories(jpegultrascan_core
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(jpegultrascan_core PUBLIC JPEG::JPEG Threads::Threads)
# public: code that includes jpeglib.h next to ours needs to know which libjpeg it gets
if(JPEGULTRASCAN_IS_MOZJPEG)
    target_compile_definitions(jpegultrascan_core PUBLIC JPEGULTRASCAN_MOZJPEG=1)
endif()
# chisel builds its release in unity mode; these sources gain nothing from it
set_target_properties(jpegultrascan_core PROPERTIES UNITY_BUILD OFF)
add_library(jpegultrascan::core ALIAS jpegultrascan_core)

if(JPEGULTRASCAN_BUILD_CLI)
    add_executable(jpegultrascancpp src/main.cpp)
    target_link_libraries(jpegultrascancpp PRIVATE jpegultrascan_core)
    target_include_directories(jpegultrascancpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_compile_definitions(jpegultrascancpp PRIVATE JPEGULTRASCAN_VERSION="${PROJECT_VERSION}")
endif()

# the tests aren't part of the repository; a checkout that has them builds them
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt")
    add_subdirectory(tests)
endif()
