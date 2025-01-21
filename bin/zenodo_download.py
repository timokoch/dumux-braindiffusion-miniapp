# SPDX-FileCopyrightText: Copyright © Jørgen Riseth and Timo Koch
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Download data from Zenodo.
In case the data is not publicly avaiable yet, set the
access token in the environment e.g. by calling:
ZENODO_ACCESS_TOKEN=<your_token> python zenodo_download.py
"""


import argparse
import requests
from pathlib import Path
import os


class Settings:
    def __init__(self):
        self.api_url = "https://zenodo.org/api/deposit/depositions/14266867"
        self.access_token = os.environ.get("ZENODO_ACCESS_TOKEN", None)
        if self.access_token:
            print("-- Found access token in environment variable ZENODO_ACCESS_TOKEN")


def get_file_list(settings: Settings) -> list[dict[str, str]]:
    r = requests.get(
        f"{settings.api_url}/files",
        params={"access_token": settings.access_token} if settings.access_token else {},
    )
    r.raise_for_status()
    return r.json()


def find_file_id(filename: str, files: list[dict[str, str]], settings: Settings) -> str:
    file_ids = [
        deposition_file["id"]
        for deposition_file in files
        if deposition_file["filename"] == filename
    ]
    if len(file_ids) == 0:
        raise ValueError(f"Couldn't find {filename} in deposition")
    if len(file_ids) > 1:  # Not sure if this can happen.
        raise ValueError(f"Multiple files found named {filename}.")
    return file_ids[0]


def download_single_file(file_download_url: str, output: str, settings: Settings):
    r = requests.get(
        file_download_url,
        params={"access_token": settings.access_token} if settings.access_token else {},
    )
    r.raise_for_status()
    with open(output, "wb") as f:
        f.write(r.content)


def list_all_files(files: list[dict[str, str]]):
    for deposition_file in files:
        print(deposition_file["filename"])


def validate_input(parser: argparse.ArgumentParser):
    args = parser.parse_args()
    if args.ls:
        return args
    elif not args.output:
        raise parser.error("Missing required option '--output'.")
    if (not args.all) and (args.filename is None):
        parser.error(
            "Missing '--filename'. Option required without '--list' or '--all'"
        )
    return args


def main():
    parser = argparse.ArgumentParser(description="Download files.")
    parser.add_argument(
        "--output",
        type=str,
        help="Output directory.",
    )
    parser.add_argument(
        "--filename",
        type=str,
        help="Name of file to download. Ignored if '--all' is used",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Download all available files.",
    )
    parser.add_argument(
        "--list",
        dest="ls",
        action="store_true",
        help="List files in record, without downloading.",
    )

    args = validate_input(parser)

    settings = Settings()
    files = get_file_list(settings)

    if args.ls:
        list_all_files(files)
        return

    if not args.all:
        if args.filename is None:
            parser.error("'--filename' option required without '--list' or '--all'")
        files = [file for file in files if file["filename"] == args.filename]
        if len(files) == 0:
            raise ValueError(f"Couldn't find {args.filename} in deposition")

    output_dir = Path(args.output)
    if not output_dir.exists():
        raise IOError(f"Output directory {output_dir} does not exist")

    for file in files:
        output_name = output_dir / Path(file["filename"])
        if output_name.is_file():
            print(f"-- Found file {file['filename']}. Skip download.")
            continue

        print(f"Downloading file {file['filename']}")
        download_single_file(
            file["links"]["download"], f"{args.output}/{file['filename']}", settings
        )


if __name__ == "__main__":
    main()
