# SPDX-FileCopyrightText: Copyright © Timo Koch
# SPDX-License-Identifier: CC0-1.0

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

find_package(Python COMPONENTS Interpreter REQUIRED)
message(STATUS "Found Python interpreter: ${Python_EXECUTABLE}")

set(MRI_CONCENTRATION_DATA_FILE ${CMAKE_CURRENT_SOURCE_DIR}/data/mesh-data/mri_processed_data/sub-01/modeling/resolution32/data.vtk)
set(MRI_T1_DATA_FILE ${CMAKE_CURRENT_SOURCE_DIR}/data/mri-dataset/mri_dataset/sub-01/ses-01/anat/sub-01_ses-01_T1w.nii.gz)

add_custom_command(
  OUTPUT ${MRI_CONCENTRATION_DATA_FILE} ${MRI_T1_DATA_FILE}
  COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/bin/get_data.py --output ${CMAKE_CURRENT_SOURCE_DIR}/data
  BYPRODUCTS ${CMAKE_CURRENT_SOURCE_DIR}/data/mesh-data.zip
             ${CMAKE_CURRENT_SOURCE_DIR}/data/mri-dataset.zip
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_custom_target(
  download_mri_data ALL
  DEPENDS ${MRI_CONCENTRATION_DATA_FILE} ${MRI_T1_DATA_FILE}
)

# the application code is in the app/ subdirectory
add_subdirectory(app)
add_subdirectory(doc)
add_subdirectory(cmake/modules)

finalize_dune_project()
