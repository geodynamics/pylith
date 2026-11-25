#!/usr/bin/env python3
# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


"""Fix the preamble and the `End of file` marker in source files.

The fixes mirror the checks in etc/ci/conventions/check.py:

  * Every file opens with the standard preamble (mode line, coding line, project text,
    copyright line, and SPDX license line), written with the comment syntax of the language
    of the file. A missing preamble is added, and a preamble with the standard text but a wrong
    layout is rewritten, using PREAMBLE_TEMPLATE. A preamble with other text (such as a
    legacy copyright notice) is reported and left unchanged.
  * The preamble is followed by the number of blank lines the language calls for.
  * The file closes with a single `End of file` marker, preceded by the number of blank lines
    the language calls for, and followed by exactly one newline.

By default, files are fixed in place. Use --dry-run to only list files that need changes.

Example:
    fix-preambles-tails.py --path libsrc --regex '.*\\.(cc|hh|icc)$' --dry-run

A missing marker is added. Files with more than one marker, or with text after the marker,
cannot have their marker fixed automatically. They are reported as errors; any other fix to
them is still applied.

Files matching the patterns in the `exempt` list of [tool.conventions] in pyproject.toml are
skipped, like check.py skips them. The file is looked for in --path and the directories above it.

Exit status is 0 if no file needs changes (or all were fixed) and 1 if any file needs changes in
--dry-run mode or could not be fixed.
"""

import argparse
import dataclasses
import datetime
import fnmatch
import pathlib
import re
import sys
import tomllib
from collections.abc import Iterator, Sequence

# The standard preamble. Each line is formatted with {c} (the comment character of the
# language), {mode} (the editor mode of the language), {year} (the last year of the copyright),
# and {LICENSE_TAG} (the SPDX identifier of the license). A line holding only {c} is a spacer.
PREAMBLE_TEMPLATE = """\
{c} -*- {mode} -*-
{c} -*- coding: utf-8 -*-
{c}
{c} This code is part of PyLith, developed through the Computational Infrastructure
{c} for Geodynamics (https://github.com/geodynamics/pylith).
{c}
{c} Copyright (c) 2010-{year} The PyLith authors and contributors
{c} SPDX-License-Identifier: {LICENSE_TAG}
"""

# The license written in a new preamble.
DEFAULT_LICENSE = "MIT"
# The text of the end of file marker, after the comment character.
MARKER_TEXT = "End of file"
# The number of lines, from the top of a file, searched for the copyright line.
WINDOW = 15
# The name of the file that lists the exempt patterns, in the `exempt` list of [tool.conventions].
CONFIG_NAME = "pyproject.toml"
# The license line, which ends a preamble, with the SPDX identifier of the license.
LICENSE_RE = re.compile(r"SPDX-License-Identifier:(?:\s*(?P<tag>\S.*?))?\s*$")
# The longest comment block that is taken to be a preamble, so that a long comment that happens
# to follow the copyright line is never swallowed.
MAX_PREAMBLE = 20
# The copyright line, with any text before it and any closing year.
COPYRIGHT_RE = re.compile(r"Copyright \(c\) 2010-(?P<year>\d{4})\b")
# Any copyright line, including the legacy forms that do not match the standard text.
LEGACY_COPYRIGHT_RE = re.compile(r"Copyright \(c\)")
# The mode line a stylesheet keeps above its comment block.
OUTER_MODE_RE = re.compile(r"^/\*\s+-\*- .+ -\*-\s+\*/$")
# The text that marks an editor variable line, such as the mode and coding lines.
EDITOR_VARIABLES = "-*-"
# The lines that stay above the preamble in a file with a shebang line.
SHEBANG = r"#!"
# The lines that stay above the preamble in a web page.
WEB_PRELUDE = r"(?i:<!doctype |<\?xml )"
# The prefix of the line a c++ header uses to keep its code guard next to the preamble.
CODE_GUARD = "code guard"


class FixError(Exception):
    """A file that cannot be fixed automatically, with the reason in the message."""


@dataclasses.dataclass(frozen=True)
class Language:
    """A language of source file, as it is written in a preamble.

    The language {mode} is the editor mode named in the preamble and {comment} is the text that
    starts each preamble line. The file is recognized by one of its {suffixes}, {names}, or the
    {directories} that hold it. The {blank_after_preamble} lines follow the preamble, and one of
    the {blank_before_marker} counts of blank lines precedes the end of file marker. The lines
    that stay above the preamble match {prelude}. A language with a {block} comment syntax, a pair
    of opening and closing delimiters, puts the preamble in a block whose lines start with the
    {comment}. A stylesheet keeps its mode line outside the block, in the {outer_mode} format.
    """

    mode: str
    comment: str
    suffixes: tuple[str, ...] = ()
    names: tuple[str, ...] = ()
    directories: tuple[str, ...] = ()
    blank_after_preamble: int = 2
    blank_before_marker: tuple[int, ...] = (2,)
    block: tuple[str, str] | None = None
    outer_mode: str | None = None
    prelude: str = SHEBANG

    @property
    def marker(self) -> str:
        """Return the end of file marker of the language, with its comment syntax."""
        if self.block:
            return f"{self.block[0]} {MARKER_TEXT} {self.block[1]}"
        return f"{self.comment} {MARKER_TEXT}"

    def matches(self, path: pathlib.Path) -> bool:
        """Return whether the file at {path} is written in this language."""
        return path.name in self.names or path.suffix in self.suffixes or path.parent.name in self.directories

    def is_comment(self, line: str) -> bool:
        """Return whether the {line} is a line comment, or a spacer holding only the comment character."""
        return line.rstrip() == self.comment or line.startswith(self.comment + " ")


# The languages, in the order they are tried. Black keeps one blank line after imports, so a
# python file may close with one or two; the formatters of markdown, yaml, and python stubs
# keep a single blank line. Dockerfiles are recognized last, since any file in a docker
# directory is taken for one if it is not another language.
LANGUAGES = (
    Language("c++", "//", suffixes=(".h", ".hh", ".icc", ".cc", ".cpp", ".cu", ".i")),
    Language("c", " *", suffixes=(".c",), block=("/*", "*/")),
    Language("css", " *", suffixes=(".css",), block=("/*", "*/"), outer_mode="/*  -*- {mode} -*-  */"),
    Language("html", "", suffixes=(".html", ".htm"), block=("<!--", "-->"), prelude=WEB_PRELUDE),
    Language("python", "#", suffixes=(".py",), blank_before_marker=(2, 1)),
    Language("python", "#", suffixes=(".pyi",), blank_after_preamble=1, blank_before_marker=(1,)),
    Language("makefile", "#", suffixes=(".am",), names=("Makefile",)),
    Language("mesh", "//", suffixes=(".mesh",)),
    Language("config", "#", suffixes=(".cfg",)),
    Language("journal", "#", suffixes=(".jou",)),
    Language(
        "yaml",
        "#",
        suffixes=(".yml", ".yaml"),
        names=(".clang-format",),
        blank_after_preamble=1,
        blank_before_marker=(1,),
    ),
    Language("markdown", "%", suffixes=(".md",), blank_after_preamble=1, blank_before_marker=(1,)),
    Language("tex", "%", suffixes=(".tex",)),
    Language("autoconf", "dnl", suffixes=(".ac", ".m4")),
    Language("sh", "#", suffixes=(".sh",)),
    Language("dockerfile", "#", suffixes=(".dockerfile",), names=("Dockerfile",), directories=("docker",)),
)


def find_language(path: pathlib.Path) -> Language | None:
    """Return the language of the file at {path}, or None if it is not one of the LANGUAGES."""
    return next((language for language in LANGUAGES if language.matches(path)), None)


@dataclasses.dataclass(frozen=True)
class PreambleSpan:
    """The place of an existing preamble in the lines of a file.

    The preamble occupies the lines from {start} up to, but not including, {end}. It records
    the closing {year} of its copyright and the SPDX identifier of its {license}, either of
    which may be None.
    """

    start: int
    end: int
    year: str | None
    license: str | None = None


class PreambleTemplate:
    """The standard preamble, as text that is formatted for a language and a year."""

    def __init__(self, text: str = PREAMBLE_TEMPLATE) -> None:
        """Create the template from the {text}, one preamble line per line of text."""
        self.text = text

    def render(self, language: Language, year: str, license_tag: str = DEFAULT_LICENSE) -> list[str]:
        """Return the lines of the preamble for the {language}, with a copyright ending in {year} and the {license_tag}.

        The comment syntax and editor mode are those of the language. Trailing space is
        removed, so spacer lines hold only the comment character. A block language has its
        delimiters around the lines, and a stylesheet its mode line above the block.
        """
        values = {"c": language.comment, "mode": language.mode, "year": year, "LICENSE_TAG": license_tag}
        lines = []
        outer = []
        for template in self.text.splitlines():
            line = template.format(**values).rstrip()
            # a comment of nothing starts no line, so the space after it goes too
            line = line[1:] if not language.comment and line.startswith(" ") else line
            if language.outer_mode and "{mode}" in template:
                outer.append(language.outer_mode.format(**values))
            else:
                lines.append(line)
        if language.block:
            # a closer lines up with the stars of a comment that has them
            indent = " " if language.comment.startswith(" ") else ""
            lines = [language.block[0], *lines, indent + language.block[1]]
        return outer + lines


class PreambleLocator:
    """The search for an existing preamble in the lines of a file."""

    def __init__(self, language: Language) -> None:
        """Create the locator for files written in the {language}."""
        self.language = language
        self.prelude_re = re.compile(language.prelude)

    def prelude(self, lines: Sequence[str]) -> int:
        """Return the number of lines at the top of the {lines} that stay above the preamble."""
        count = 0
        while count < len(lines) and self.prelude_re.match(lines[count]):
            count += 1
        return count

    def locate(self, lines: Sequence[str]) -> PreambleSpan | None:
        """Return the span of the preamble in the {lines}, or None if the file has none.

        A line comment preamble is the run of comment lines around the first copyright line in
        the top of the file, below the lines that stay above it. A block comment preamble is
        the block around the copyright line, with the mode line of a stylesheet above it. A
        copyright line in any form counts, so a legacy preamble is found and not duplicated. A
        preamble that is too long, or a block that cannot be delimited, raises FixError.
        """
        top = self.prelude(lines)
        found = next(
            (
                index
                for index in range(top, min(len(lines), WINDOW))
                if LEGACY_COPYRIGHT_RE.search(lines[index])
                and (self.language.block or self.language.is_comment(lines[index]))
            ),
            None,
        )
        if found is None:
            return None
        start, end = self._block(lines, found, top) if self.language.block else self._run(lines, found, top)
        if end - start > MAX_PREAMBLE:
            raise FixError(f"comment of {end - start} lines around the copyright is too long to be a preamble")
        year = next((match["year"] for line in lines[start:end] if (match := COPYRIGHT_RE.search(line))), None)
        license = next((match["tag"] for line in lines[start:end] if (match := LICENSE_RE.search(line))), None)
        return PreambleSpan(start, end, year, license)

    def _run(self, lines: Sequence[str], found: int, top: int) -> tuple[int, int]:
        """Return the start and end of the run of line comments in the {lines} around the line {found}, below {top}."""
        start = found
        while start > top and self.language.is_comment(lines[start - 1]):
            start -= 1
        end = found + 1
        while end < len(lines) and self.language.is_comment(lines[end]):
            end += 1
        # the preamble ends with the license, and any spacers after it, and not with the comments that follow
        license = next((i for i in range(found, end) if LICENSE_RE.search(lines[i])), None)
        if license is not None:
            stop = license + 1
            while stop < end and lines[stop].rstrip() == self.language.comment:
                stop += 1
            end = stop
        return start, end

    def _block(self, lines: Sequence[str], found: int, top: int) -> tuple[int, int]:
        """Return the start and end of the block comment in the {lines} around the line {found}, below {top}.

        The delimiters must each be on a line of their own, or FixError is raised. The mode
        line of a stylesheet is part of the block when it is on the line above.
        """
        opener, closer = self.language.block
        start = next((i for i in range(found - 1, top - 1, -1) if lines[i].strip() == opener), None)
        end = next((i for i in range(found + 1, len(lines)) if lines[i].strip() == closer), None)
        if start is None or end is None:
            raise FixError(f"copyright is not in a block comment with '{opener}' and '{closer}' on lines of their own")
        if self.language.outer_mode and start > top and OUTER_MODE_RE.match(lines[start - 1]):
            start -= 1
        return start, end + 1


@dataclasses.dataclass
class Fix:
    """The outcome of fixing the text of one file.

    The fixed {lines} replace the original, {changes} describes what was changed, and the
    {errors} describe what could not be fixed.
    """

    lines: list[str]
    changes: list[str] = dataclasses.field(default_factory=list)
    errors: list[str] = dataclasses.field(default_factory=list)


class SourceFixer:
    """The fixes for the files written in one language."""

    def __init__(
        self,
        language: Language,
        template: PreambleTemplate,
        year: str,
        update_year: bool = False,
    ) -> None:
        """Create the fixer for the {language}, writing preambles from the {template}.

        A new preamble ends its copyright in {year}. The year of an existing preamble is kept,
        unless {update_year} is set.
        """
        self.language = language
        self.template = template
        self.year = year
        self.update_year = update_year
        self.locator = PreambleLocator(language)
        self.marker_re = re.compile(r"^[ \t]*" + re.escape(language.marker) + r"[ \t]*$")

    def fix_preamble(self, lines: list[str]) -> Fix:
        """Return the {lines} with the standard preamble, and the blank lines that follow it.

        A missing preamble is added below the lines that stay above it. A preamble with the
        standard text is rewritten when it differs from the template in anything but trailing
        space; the editor variable lines (mode and coding) are not part of the text, and the
        license, whatever it is, is kept. A preamble with any other text is kept, reported as an error, and only its blank lines are fixed.
        """
        span = self.locator.locate(lines)
        if span is None:
            start = end = self.locator.prelude(lines)
            year, license = self.year, DEFAULT_LICENSE
        else:
            start, end = span.start, span.end
            year = self.year if self.update_year or not span.year else span.year
            license = span.license or DEFAULT_LICENSE
        preamble = self.template.render(self.language, year, license)
        body = lines[end:]
        while body and not body[0].strip():
            body.pop(0)
        blanks = self.language.blank_after_preamble
        if body and body[0].startswith(f"{self.language.comment} {CODE_GUARD}"):
            blanks = 1
        head = lines[:start]
        # a block sits a blank line below the lines above it
        if span is None and head and self.language.block:
            head = head + [""]
        current = lines[start:end]
        fix = Fix([])
        if span is None:
            middle = preamble
            fix.changes.append("added preamble")
        else:
            # the text is standard, whatever the year of the copyright and the editor variables
            shown = [line.rstrip() for line in current]
            # spacers after the license are not part of the text
            while shown and shown[-1] == self.language.comment and not self.language.block:
                shown.pop()
            existing = self.template.render(self.language, span.year or year, license)
            standard = self._text(shown) == self._text(existing)
            if not standard:
                middle = current
                fix.errors.append("non-standard preamble text; preamble left unchanged")
            elif shown != preamble or len(shown) != len(current):
                middle = preamble
                fix.changes.append("rewrote preamble")
            else:
                # differing only in trailing space, it is left as it is
                middle = current
        fix.lines = head + middle + [""] * blanks + body
        if not fix.changes and lines[end : end + blanks + 1] != [""] * blanks + body[:1]:
            fix.changes.append("fixed blank lines after preamble")
        return fix

    @staticmethod
    def _text(lines: Sequence[str]) -> list[str]:
        """Return the {lines} of a preamble without its editor variable lines."""
        return [line for line in lines if EDITOR_VARIABLES not in line]

    def _locate(self, lines: Sequence[str]) -> PreambleSpan | None:
        """Return the span of the preamble in the {lines}, or None if there is none or it cannot be delimited."""
        try:
            return self.locator.locate(lines)
        except FixError:
            return None

    def fix_marker(self, lines: list[str]) -> Fix:
        """Return the {lines} with the blank lines before the end of file marker that the language calls for.

        A file with no marker gets one. A file with more than one marker, or text after the
        marker, is returned with an error. The blank lines are left alone when the file has no
        body, since then they are the ones that follow the preamble, and when their count is one
        the language allows.
        """
        fix = Fix(list(lines))
        marks = [index for index, line in enumerate(lines) if self.marker_re.match(line)]
        if not marks:
            return self.add_marker(lines)
        if len(marks) > 1:
            fix.errors.append(f"multiple '{self.language.marker}' markers")
            return fix
        body = lines[: marks[0]]
        if any(line.strip() for line in lines[marks[0] + 1 :]):
            fix.errors.append(f"text found after '{self.language.marker}' marker")
            return fix
        span = self._locate(lines)
        floor = span.end if span else 0
        blanks = len(body) - len(list(self._strip(body)))
        if not any(line.strip() for line in body[floor:]):
            # no body: the blank lines belong to the preamble
            return fix
        if blanks in self.language.blank_before_marker:
            kept = lines[: marks[0]]
        else:
            kept = list(self._strip(body)) + [""] * self.language.blank_before_marker[0]
            fix.changes.append("fixed blank lines before end of file marker")
        fix.lines = kept + [self.language.marker]
        return fix

    def add_marker(self, lines: list[str]) -> Fix:
        """Return the {lines} with the end of file marker added at the end.

        The marker follows the blank lines the language calls for. They are those before the
        marker when the file has a body, and those after the preamble when it has none.
        """
        kept = list(self._strip(lines))
        span = self._locate(kept)
        has_body = any(line.strip() for line in kept[span.end if span else 0 :])
        if has_body:
            blanks = self.language.blank_before_marker[0]
        else:
            blanks = self.language.blank_after_preamble if span else 0
        fix = Fix(kept + [""] * blanks + [self.language.marker])
        fix.changes.append("added end of file marker")
        return fix

    @staticmethod
    def _strip(lines: Sequence[str]) -> Iterator[str]:
        """Yield the {lines} without the blank lines at the end."""
        end = len(lines)
        while end and not lines[end - 1].strip():
            end -= 1
        yield from lines[:end]

    def fix_text(self, text: str, preamble: bool = True, marker: bool = True) -> tuple[str, Fix]:
        """Return the fixed {text} of a file, and the account of the fixes.

        The preamble is fixed when {preamble} is set, and the end of file marker is fixed when
        {marker} is set. The line endings of the text are kept.
        """
        newline = "\r\n" if "\r\n" in text else "\n"
        lines = text.replace("\r\n", "\n").split("\n")
        if lines and lines[-1] == "":
            lines.pop()
        total = Fix(lines)
        if preamble:
            try:
                step = self.fix_preamble(total.lines)
            except FixError as err:
                total.errors.append(str(err))
            else:
                total.lines = step.lines
                total.changes += step.changes
                total.errors += step.errors
        if marker:
            step = self.fix_marker(total.lines)
            total.lines = step.lines
            total.changes += step.changes
            total.errors += step.errors
        fixed = newline.join(total.lines) + newline if total.lines else ""
        if fixed != text and not total.changes and not total.errors:
            total.changes.append("fixed trailing newlines")
        return fixed, total


def find_config(path: pathlib.Path) -> pathlib.Path | None:
    """Return the pyproject.toml in the directory {path} or the nearest directory above it, or None if there is none."""
    start = path.resolve()
    return next((d / CONFIG_NAME for d in (start, *start.parents) if (d / CONFIG_NAME).is_file()), None)


def load_exemptions(config: pathlib.Path) -> tuple[str, ...]:
    """Return the patterns of the exempt files listed in the {exempt} list of [tool.conventions] in the file {config}."""
    with config.open("rb") as stream:
        table = tomllib.load(stream)
    return tuple(table.get("tool", {}).get("conventions", {}).get("exempt", []))


@dataclasses.dataclass
class Tally:
    """The count of the files checked, changed, and in error."""

    checked: int = 0
    changed: int = 0
    errors: int = 0


class Runner:
    """The fixing of every file under a directory."""

    def __init__(
        self,
        path: pathlib.Path,
        pattern: re.Pattern[str],
        root: pathlib.Path,
        exclude: Sequence[str],
        fixers: dict[Language, SourceFixer],
        dry_run: bool,
        preamble: bool,
        marker: bool,
    ) -> None:
        """Create the runner for the files under {path} whose full path matches {pattern}.

        Files whose path relative to the repository {root} matches a glob in {exclude} are
        skipped; these are the exempt patterns of the repository and any given on the command
        line. Each file
        is fixed by the fixer of its language in {fixers}; the files are only reported when
        {dry_run} is set, and the {preamble} and {marker} fixes can each be turned off.
        """
        self.path = path
        self.pattern = pattern
        self.root = root.resolve()
        self.exclude = exclude
        self.fixers = fixers
        self.dry_run = dry_run
        self.preamble = preamble
        self.marker = marker
        self.tally = Tally()

    def find_files(self) -> Iterator[pathlib.Path]:
        """Yield the files under the path to fix: those that match the pattern, are not excluded, and are not in .git."""
        for filename in sorted(self.path.rglob("*")):
            if (
                filename.is_file()
                and ".git" not in filename.relative_to(self.path).parts
                and self.pattern.search(str(filename))
                and not self.excluded(filename)
            ):
                yield filename

    def excluded(self, filename: pathlib.Path) -> bool:
        """Return whether the file at {filename} matches an exclude pattern, by its path relative to the root."""
        try:
            relative = filename.resolve().relative_to(self.root)
        except ValueError:
            relative = filename.relative_to(self.path)
        return any(fnmatch.fnmatch(relative.as_posix(), glob) for glob in self.exclude)

    def process(self, filename: pathlib.Path) -> None:
        """Fix the file at {filename}, report what was done, and update the tally.

        Files in a language that is not known, and files that are not text, are skipped.
        """
        language = find_language(filename)
        if language is None:
            return
        try:
            # newline="" preserves existing line endings.
            with filename.open("r", newline="", encoding="utf-8") as f:
                text = f.read()
        except UnicodeDecodeError:
            return
        self.tally.checked += 1
        new_text, fix = self.fixers[language].fix_text(text, self.preamble, self.marker)
        for error in fix.errors:
            print(f"{filename}: ERROR: {error}")
        if fix.errors:
            self.tally.errors += 1
        if new_text == text:
            return
        self.tally.changed += 1
        details = ", ".join(fix.changes)
        if self.dry_run:
            print(f"{filename}: needs fix ({details})")
        else:
            with filename.open("w", newline="", encoding="utf-8") as f:
                f.write(new_text)
            print(f"{filename}: fixed ({details})")

    def run(self) -> int:
        """Fix every file, print the summary, and return the exit status: 1 if any file needs changes or is in error."""
        for filename in self.find_files():
            self.process(filename)
        verb = "need changes" if self.dry_run else "fixed"
        print(f"Checked {self.tally.checked} file(s): {self.tally.changed} {verb}, {self.tally.errors} error(s).")
        return 1 if (self.tally.errors or (self.dry_run and self.tally.changed)) else 0


def build_parser() -> argparse.ArgumentParser:
    """Return the parser of the command line arguments."""
    parser = argparse.ArgumentParser(
        description=f"Fix the preamble and the '{MARKER_TEXT}' marker of source files.",
    )
    parser.add_argument(
        "--regex",
        default=".",
        help="Regular expression matched (re.search) against file paths. Default: %(default)s",
    )
    parser.add_argument(
        "--path",
        default=pathlib.Path("."),
        type=pathlib.Path,
        help="Directory to search recursively. Default: %(default)s",
    )
    parser.add_argument(
        "--exclude",
        action="append",
        default=[],
        metavar="GLOB",
        help="Glob matched against the path relative to the repository root; matching files are skipped, "
        "in addition to the exempt patterns in pyproject.toml. May be repeated.",
    )
    parser.add_argument(
        "--config",
        type=pathlib.Path,
        help=f"The {CONFIG_NAME} with the exempt patterns; its directory is the repository root. "
        "Default: the first one in --path or the directories above it.",
    )
    parser.add_argument(
        "--only",
        choices=("preamble", "marker"),
        help="Fix only the preamble or only the end of file marker. Default: both.",
    )
    parser.add_argument(
        "--year",
        default=str(datetime.date.today().year),
        help="Last year of the copyright in new preambles. Default: %(default)s",
    )
    parser.add_argument(
        "--update-year",
        action="store_true",
        help="Also set the year of existing preambles to --year.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only show files needing changes; do not modify them.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Fix the files named by the command line {argv}, or by sys.argv, and return the exit status."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        pattern = re.compile(args.regex)
    except re.error as err:
        parser.error(f"invalid regular expression: {err}")
    if not re.fullmatch(r"\d{4}", args.year):
        parser.error(f"year must have four digits: {args.year}")
    if not args.path.is_dir():
        parser.error(f"search path is not a directory: {args.path}")

    config = args.config or find_config(args.path)
    if config is not None and not config.is_file():
        parser.error(f"configuration file not found: {config}")
    exempt = load_exemptions(config) if config else ()
    root = config.parent if config else args.path

    template = PreambleTemplate()
    fixers = {language: SourceFixer(language, template, args.year, args.update_year) for language in LANGUAGES}
    runner = Runner(
        path=args.path,
        pattern=pattern,
        root=root,
        exclude=[*exempt, *args.exclude],
        fixers=fixers,
        dry_run=args.dry_run,
        preamble=args.only != "marker",
        marker=args.only != "preamble",
    )
    return runner.run()


if __name__ == "__main__":
    sys.exit(main())


# End of file
