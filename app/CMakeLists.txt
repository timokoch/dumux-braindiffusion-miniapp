# SPDX-FileCopyrightText: Copyright © Timo Koch
# SPDX-License-Identifier: CC0-1.0

dune_symlink_to_source_files(FILES params.input curvefit.input)

dumux_add_test(
    NAME braindiffusion_test_compile
    SOURCES main.cc
    COMPILE_ONLY
)

dumux_add_test(
    NAME braindiffusion
    SOURCES main.cc
)

add_dependencies(braindiffusion download_mri_data)
dune_target_enable_all_packages(braindiffusion)

dumux_add_test(
    NAME curvefit
    SOURCES curvefit.cc
)

add_dependencies(curvefit download_mri_data)
dune_target_enable_all_packages(curvefit)
