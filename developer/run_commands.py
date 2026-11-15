#!/usr/bin/env python3
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
"""Run commands listed in a (YAML) file.

This script is primarily used to run an example suite or a benchmark, including pre- and post-processing scripts.
"""

import pathlib
import subprocess

import yaml


class Runner:

    def run(self, filename: str, dir_path: pathlib.Path, log_path: pathlib.Path) -> int:
        self._run_dir(filename=filename, dir_path=dir_path, log_path=log_path)
        for _, dir_names, _ in dir_path.walk():
            for subdir_name in sorted(dir_names):
                self._run_dir(filename=filename, dir_path=pathlib.Path(subdir_name), log_path=log_path)
        return 0

    def _get_commands(self, file_path: pathlib.Path) -> list[str]:
        """Get commands from YAML file {file_path} to run example suite."""
        with open(file_path, encoding="utf-8") as fin:
            config = yaml.safe_load(fin)
        if "commands" in config:
            return config["commands"]
        else:
            print(f"Skipping {file_path}. Missing 'commands'.")
        return []

    def _run_commands(self, subdir: pathlib.Path, log_path: pathlib.Path, commands: list[str]) -> None:
        """Run {commands} in {subdir}."""
        for cmd in commands:
            result = subprocess.run(cmd, cwd=subdir, shell=True, capture_output=True, text=True, check=False)
            if result.returncode != 0:
                print(f"ERROR ({result.returncode}): {subdir} {cmd}")
                print(result.stderr.strip())
            else:
                print(f"OK: {subdir} {cmd}")
                if log_path:
                    with open(log_path, mode="a", encoding="utf-8") as output:
                        output.write(f"\nRUNNING {subdir} {cmd}\n")
                        output.write(result.stdout)

    def _run_dir(self, filename: str, dir_path: pathlib.Path, log_path: pathlib.Path) -> None:
        """Run example suite under {path}."""
        config_path = dir_path / filename
        if config_path.is_file():
            commands = self._get_commands(config_path)
            if commands:
                self._run_commands(subdir=dir_path, log_path=log_path, commands=commands)


def cli() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--path", dest="path", default=".", help="Path to search for {filename}.")
    parser.add_argument("--filename", dest="filename", default="examples.yaml", help="YAML file with list of commands.")
    parser.add_argument("--log", dest="log", help="Write stdout to {log} file.")
    args = parser.parse_args()

    dir_path = pathlib.Path(args.path)
    log_path = None
    if args.log:
        log_path = pathlib.Path(args.log)
        log_path.unlink(missing_ok=True)
    return Runner().run(filename=args.filename, dir_path=dir_path, log_path=log_path)


if __name__ == "__main__":
    raise SystemExit(cli())


# End of file
