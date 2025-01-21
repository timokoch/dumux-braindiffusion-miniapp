# SPDX-FileCopyrightTest: Copyright © Timo Koch
# SPDX-License-Identifier: GPL-3.0-or-later

from zenodo_download import download_single_file, get_file_list, Settings
from pathlib import Path
import argparse
import shutil


def main():
    parser = argparse.ArgumentParser(description="Download files.")
    parser.add_argument(
        "--output",
        type=str,
        help="Output directory.",
    )

    args = parser.parse_args()

    settings = Settings()
    files = get_file_list(settings)
    output_dir = Path(args.output)
    if not output_dir.exists():
        raise IOError(f"Output directory {output_dir} does not exist")

    for file in files:
        file_name = file["filename"]
        if file_name not in ["mesh-data.zip", "mri-dataset.zip"]:
            continue

        output_name = output_dir / Path(file_name)
        if output_name.is_file():
            print(f"-- Found file {file_name}. Skip download.")
            continue

        print(f"-- Downloading file {file_name}")
        download_single_file(file["links"]["download"], output_name, settings)

    for file in files:
        file_name = file["filename"]
        if file_name not in ["mesh-data.zip", "mri-dataset.zip"]:
            continue

        output_unpack_dir = output_dir / Path(file_name.rstrip(".zip"))
        if output_unpack_dir.is_dir():
            print(
                f"-- Found extracted data directory {output_unpack_dir}. Skip unpacking."
            )
            continue

        output_name = output_dir / Path(file_name)
        print(f"-- Unpacking archive to {output_unpack_dir}")
        shutil.unpack_archive(output_name, output_unpack_dir)


if __name__ == "__main__":
    main()
