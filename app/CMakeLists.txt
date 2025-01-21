# SPDX-FileCopyrightText: Copyright © Timo Koch
# SPDX-License-Identifier: CC0-1.0

dune_symlink_to_source_files(FILES params.input curvefit.input data.vtu plot.py)

dumux_add_test(
    NAME braindiffusion
    SOURCES main.cc
)

add_dependencies(braindiffusion download_mri_data)

dumux_add_test(
    NAME curvefit
    SOURCES curvefit.cc
)
