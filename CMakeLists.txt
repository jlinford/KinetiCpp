cmake_minimum_required(VERSION 3.21)

project(
  KinetiCpp
  VERSION 0.0.1
  DESCRIPTION ""
  HOMEPAGE_URL ""
  LANGUAGES CXX)

# Suppress warnings about unused C compiler
set (ignoreCC ${CMAKE_C_COMPILER})

# Prevent in-source builds
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "In-source builds not allowed. Please run CMake in a build directory.")
endif()

# Are we building for ourselves, or others?
set(KINETICPP_STANDALONE_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(KINETICPP_STANDALONE_PROJECT ON)
endif()

if (KINETICPP_STANDALONE_PROJECT)

    set(CMAKE_CXX_STANDARD 23)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)

    # generate compile_commands.json
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        find_program(CLANG_TIDY_PROGRAM NAMES clang-tidy)
        if(CLANG_TIDY_PROGRAM)
            set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_PROGRAM}")
        endif()
    endif()

    find_program(CCACHE_PROGRAM ccache)
    if(CCACHE_PROGRAM)
        SET(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    endif(CCACHE_PROGRAM)

    include(cmake/CheckCxxCompilerFeatures.cmake)

    add_library(KinetiCpp INTERFACE)
    add_library(kineticpp ALIAS KinetiCpp)
    add_library(kineticpp::kineticpp ALIAS KinetiCpp)

    target_compile_features(KinetiCpp INTERFACE cxx_std_23)

    target_include_directories(KinetiCpp
      INTERFACE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include>
    )

    include(cmake/Dependencies.cmake)

    add_subdirectory(examples)

else()
  add_library(KinetiCpp INTERFACE)
  add_library(kineticpp ALIAS KinetiCpp)
  add_library(kineticpp::kineticpp ALIAS KinetiCpp)
  set_property(TARGET KinetiCpp PROPERTY CXX_STANDARD 23)
  target_include_directories(KinetiCpp PUBLIC ${PROJECT_SOURCE_DIR}/src)
endif()
