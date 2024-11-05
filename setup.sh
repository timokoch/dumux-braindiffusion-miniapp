#!/usr/bin/env bash

git clone -b releases/2.10 https://gitlab.dune-project.org/core/dune-common.git
git clone -b releases/2.10 https://gitlab.dune-project.org/core/dune-geometry.git
git clone -b releases/2.10 https://gitlab.dune-project.org/core/dune-grid.git
git clone -b releases/2.10 https://gitlab.dune-project.org/core/dune-istl.git
git clone -b releases/2.10 https://gitlab.dune-project.org/core/dune-localfunctions.git
git clone -b releases/2.10 https://gitlab.dune-project.org/extensions/dune-alugrid.git
git clone -b 3f1db4cb2c085433cc82bc68b062146285bc3c8e https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git

# configure & build libs
./dune-common/bin/dunecontrol --opts=dumux-braindiffusion-miniapp/cmake.opts cmake
./dune-common/bin/dunecontrol --opts=dumux-braindiffusion-miniapp/cmake.opts make -j
