#!/usr/bin/env python3
"""Verify that a gzipped tarball from `make dist` contains all tracked files in the repository.
 
The repository file list comes from `git ls-files` run in the current
directory, so an unclean working tree (extra/untracked files) does not affect
the result.
"""

import sys
import fnmatch
import tarfile
import argparse
import subprocess
from pathlib import Path


def repo_files():
    """Return the set of tracked repo file paths via `git ls-files`.
 
    Raises subprocess.CalledProcessError if `git ls-files` fails.
    """
    output = subprocess.run(
        ["git", "ls-files", "-z"],
        capture_output=True,
        check=True,
        text=True,
    ).stdout
    return {p for p in output.split("\0") if p}


def tarball_files(path, strip_prefix):
    """Return the set of file paths in the tarball with the prefix stripped."""
    files = set()
    with tarfile.open(path, "r:gz") as tf:
        for member in tf.getmembers():
            if not member.isfile():
                continue
            name = member.name
            if strip_prefix and name.startswith(strip_prefix):
                name = name[len(strip_prefix):]
            files.add(name)
    return files


def load_excludes(path):
    """Load exclusion patterns, ignoring blank lines and comments."""
    patterns = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            patterns.append(line)
    return patterns


def is_excluded(path, patterns):
    return any(fnmatch.fnmatch(path, pat) for pat in patterns)


def derive_strip_prefix(tarball):
    """Derive the strip prefix from a tarball filename.
 
    e.g. '/path/to/repo-x.y.z.tar.gz' -> 'repo-x.y.z'
    """
    name = Path(tarball).name
    for suffix in (".tar.gz", ".tgz"):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return Path(name).stem


def verify(tarball, exclude_file, strip_prefix):
    """Run the verification and return an exit code."""
    repo = repo_files()
    tar = tarball_files(tarball, strip_prefix)
    excludes = load_excludes(exclude_file)

    missing = sorted(
        p for p in repo
        if p not in tar and not is_excluded(p, excludes)
    )

    if missing:
        print("Files missing from the tarball:")
        for p in missing:
            print("  " + p)
        return 1

    print("The tarball is OK.")
    return 0


def cli(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tarball", required=True, help="Path to the gzipped tarball")
    parser.add_argument("--exclude-file", default=Path("developer") / "tarball_exclude.txt", help="File listing paths/patterns to exclude")
    parser.add_argument("--strip-prefix", default=None,
        help="leading path prefix to strip from tarball member names "
             "(e.g. 'repo-1.0.0'). If omitted, derived from the tarball "
             "filename (e.g. 'repo-x.y.z.tar.gz' -> 'repo-x.y.z').",
    )
    args = parser.parse_args(argv)

    strip_prefix = args.strip_prefix
    tarball_filepath = Path(args.tarball).expanduser()
    if strip_prefix is None:
        strip_prefix = derive_strip_prefix(args.tarball)
    if strip_prefix and not strip_prefix.endswith("/"):
        strip_prefix += "/"

    return verify(tarball_filepath, Path(args.exclude_file), strip_prefix)


if __name__ == "__main__":
    sys.exit(cli())
