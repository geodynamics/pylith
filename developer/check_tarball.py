#!/usr/bin/env python3
"""Verify that a gzipped tarball from `make dist` contains all tracked files in the repository.

The repository file list comes from `git ls-files` run in the current
directory (when --use-git is set), so an unclean working tree (extra/untracked
files) does not affect the result. Without --use-git, all files found under
the current directory are used instead (useful when the .git directory is not
available).
"""

from __future__ import annotations  # Needed for Python <3.10

import sys
import fnmatch
import tarfile
import argparse
import subprocess
from pathlib import Path


def repo_files() -> set[str]:
    """Return the set of tracked repo file paths via `git ls-files`.

    Paths are returned as reported by git, which uses forward slashes on all
    platforms.

    Raises:
        subprocess.CalledProcessError: If `git ls-files` fails (e.g. no .git
            directory is present).
    """
    output = subprocess.run(
        ["git", "ls-files", "-z"],
        capture_output=True,
        check=True,
        text=True,
    ).stdout
    return {p for p in output.split("\0") if p}


def local_files() -> set[str]:
    """Return the set of file paths under the current directory.

    Paths are returned relative to the current directory (equivalent to the
    format produced by `git ls-files`), using forward slashes on all
    platforms. Files inside hidden directories (e.g. .git) are skipped, but
    hidden files themselves (e.g. .gitignore) are included.
    """
    cwd = Path(".")
    files = set()
    for p in cwd.rglob("*"):
        # Skip files inside hidden directories, but allow hidden files.
        if any(part.startswith(".") for part in p.parent.parts):
            continue
        if p.is_file():
            files.add(p.as_posix())
    return files


def tarball_files(path: Path, strip_prefix: str) -> set[str]:
    """Return the set of file paths in the tarball with the prefix stripped.

    Args:
        path: Path to the gzipped tarball.
        strip_prefix: Leading path prefix to strip from each member name
            (e.g. ``'repo-1.0.0/'``). Members that do not start with the
            prefix are included as-is.

    Returns:
        A set of file paths (directories and other non-file members are
        excluded).
    """
    files = set()
    with tarfile.open(path, "r:gz") as tf:
        for member in tf.getmembers():
            if not member.isfile():
                continue
            name = member.name
            if strip_prefix and name.startswith(strip_prefix):
                name = name[len(strip_prefix) :]
            files.add(name)
    return files


def load_excludes(path: Path) -> list[str]:
    """Load exclusion patterns from a file, ignoring blank lines and comments.

    Each non-blank, non-comment line is treated as a ``fnmatch``-style glob
    pattern. Lines beginning with ``#`` are treated as comments.

    Args:
        path: Path to the exclusion patterns file.

    Returns:
        A list of glob pattern strings.
    """
    patterns = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            patterns.append(line)
    return patterns


def is_excluded(path: str, patterns: list[str]) -> bool:
    """Return True if *path* matches any of the given fnmatch-style patterns.

    Args:
        path: The file path to test.
        patterns: A list of ``fnmatch``-style glob patterns.

    Returns:
        ``True`` if *path* matches at least one pattern, ``False`` otherwise.
    """
    return any(fnmatch.fnmatch(path, pat) for pat in patterns)


def derive_strip_prefix(tarball: str) -> str:
    """Derive the strip prefix from a tarball filename.

    Strips recognised archive suffixes (``.tar.gz``, ``.tgz``) to obtain the
    bare name, which is conventionally used as the top-level directory inside
    the tarball.

    Examples:
        ``'/path/to/repo-x.y.z.tar.gz'`` → ``'repo-x.y.z'``

    Args:
        tarball: Path to (or bare filename of) the tarball.

    Returns:
        The derived prefix string (without a trailing slash).
    """
    name = Path(tarball).name
    for suffix in (".tar.gz", ".tgz"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def verify(tarball: Path, exclude_file: Path, strip_prefix: str, use_git: bool) -> int:
    """Run the verification and return an exit code.

    Compares the set of repository files against the contents of the tarball,
    reporting any files that are missing from the tarball and not covered by
    an exclusion pattern.

    Args:
        tarball: Path to the gzipped tarball to verify.
        exclude_file: Path to the file containing exclusion patterns.
        strip_prefix: Leading prefix to strip from tarball member names.
        use_git: If ``True``, use ``git ls-files`` to enumerate repository
            files; otherwise scan the current directory with
            :func:`local_files`.

    Returns:
        ``0`` if the tarball contains all expected files, ``1`` otherwise.
    """
    repo = repo_files() if use_git else local_files()
    tar = tarball_files(tarball, strip_prefix)
    excludes = load_excludes(exclude_file)

    missing = sorted(p for p in repo if p not in tar and not is_excluded(p, excludes))

    if missing:
        print("Files missing from the tarball:")
        for p in missing:
            print("  " + p)
        return 1

    print("The tarball is OK.")
    return 0


def cli(argv: list[str] | None = None) -> int:
    """Parse command-line arguments and run the verification.

    Args:
        argv: Argument list to parse. Defaults to ``sys.argv[1:]`` when
            ``None``.

    Returns:
        Exit code: ``0`` on success, ``1`` if files are missing.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tarball", required=True, help="Path to the gzipped tarball")
    parser.add_argument(
        "--exclude-file",
        default=Path("developer") / "tarball_exclude.txt",
        help="File listing paths/patterns to exclude",
    )
    parser.add_argument(
        "--strip-prefix",
        default=None,
        help="leading path prefix to strip from tarball member names "
        "(e.g. 'repo-1.0.0'). If omitted, derived from the tarball "
        "filename (e.g. 'repo-x.y.z.tar.gz' -> 'repo-x.y.z').",
    )
    parser.add_argument(
        "--use-git",
        action="store_true",
        default=False,
        help="Use `git ls-files` to enumerate repository files instead of "
        "scanning the current directory. Requires a .git directory.",
    )
    args = parser.parse_args(argv)

    strip_prefix = args.strip_prefix
    tarball_filepath = Path(args.tarball).expanduser()
    if strip_prefix is None:
        strip_prefix = derive_strip_prefix(args.tarball)
    if strip_prefix and not strip_prefix.endswith("/"):
        strip_prefix += "/"

    return verify(tarball_filepath, Path(args.exclude_file), strip_prefix, args.use_git)


if __name__ == "__main__":
    sys.exit(cli())
