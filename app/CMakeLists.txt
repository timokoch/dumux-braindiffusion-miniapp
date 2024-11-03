dune_symlink_to_source_files(FILES params.input data.vtu plot.py)

dumux_add_test(
    NAME braindiffusion
    SOURCES main.cc
)