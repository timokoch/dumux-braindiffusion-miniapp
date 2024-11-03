cmake_minimum_required(VERSION 3.18)
project(dumux-braindiffusion-miniapp CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
  ${dune-common_MODULE_PATH})
include(DuneMacros)

dune_project()

# the application code is in the app/ subdirectory
add_subdirectory(app)
add_subdirectory(doc)
add_subdirectory(cmake/modules)

finalize_dune_project()
