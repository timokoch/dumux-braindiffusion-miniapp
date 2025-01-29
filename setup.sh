#!/usr/bin/env bash

git clone -b releases/2.10 https://gitlab.dune-project.org/core/dune-common.git
git clone -b releases/2.10 https://gitlab.dune-project.org/core/dune-geometry.git
git clone -b releases/2.10 https://gitlab.dune-project.org/core/dune-grid.git
git clone -b releases/2.10 https://gitlab.dune-project.org/core/dune-istl.git
git clone -b releases/2.10 https://gitlab.dune-project.org/core/dune-localfunctions.git
git clone -b releases/2.10 https://gitlab.dune-project.org/extensions/dune-alugrid.git
git clone -b master https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git
cd dumux && git submodule update --init && cd ..
cd dumux && git reset --hard 1e239fcb1d7a11aae5cb76d58a08b700609cfdfb && cd ..

# configure & build libs
./dune-common/bin/dunecontrol --opts=dumux-braindiffusion-miniapp/cmake.opts cmake
./dune-common/bin/dunecontrol --opts=dumux-braindiffusion-miniapp/cmake.opts make -j
