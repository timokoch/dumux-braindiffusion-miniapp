#!/usr/bin/env bash

# reconfigure & build libs
./dune-common/bin/dunecontrol --opts=dumux-braindiffusion-miniapp/cmake.opts bexec rm -r CMakeCache.txt CMakeFiles
./dune-common/bin/dunecontrol --opts=dumux-braindiffusion-miniapp/cmake.opts cmake
./dune-common/bin/dunecontrol --opts=dumux-braindiffusion-miniapp/cmake.opts make -j
