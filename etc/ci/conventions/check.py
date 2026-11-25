#!/usr/bin/env python3
# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


# externals
import argparse
import ast
import datetime
import fnmatch
import io
import os
import re
import subprocess
import sys
import tokenize
import tomllib

# the description of the program, for the help screen
DESCRIPTION = """Check a code checkout against the conventions of the repository:
  (1) spelling,
  (2) the formatting of Python and C++ sources,
  (3) file preambles and closing markers, and
  (4) commit messages,
The same checks also run on every pull request.
"""

# the checks, in the order they run
CHECKS = ("spelling", "python", "cxx", "preambles", "commits")

# black formatting options
BLACK_OPTIONS = (
    "--line-length",
    "120",
)

# the files whose preambles and closing markers are not checked, by extension, unless they carry
# a preamble anyway: binaries, and formats that cannot carry comments
OPAQUE = (
    ".png",
    ".jpg",
    ".gif",
    ".ico",
    ".svg",
    ".pdf",
    ".graffle",
    ".numbers",
    ".nb",
    ".odb",
    ".ttf",
    ".woff",
    ".woff2",
    ".dylib",
    ".json",
    ".ipynb",
    ".csv",
    ".txt",
    ".dict",
    ".dat",
    ".data",
    ".pub",
    ".plist",
)
# the c++ sources clang-format looks after
CXX = (".h", ".hh", ".icc", ".cc", ".cpp", ".cu")
# the yaml files, whose formatter keeps a single blank line wherever it finds more
YAML = (".yml", ".yaml", ".clang-format")
# the copyright line, with whatever comment marker precedes it
COPYRIGHT_LINE = re.compile(r"^(?P<prefix>.*?)Copyright \(c\) 2010-\d{4} The PyLith authors and contributors$")
# the license line, which follows the copyright; the license itself is any text
LICENSE_LINE = re.compile(r"SPDX-License-Identifier: \S.*$")
# the lines that open a block comment preamble
OPENERS = ("<!--", "{/*", "/*")
# and the ones that close it
CLOSERS = ("-->", "*/}", "*/")
# the lines that may precede the preamble: interpreters, xml declarations and doctypes, along
# with their indented continuations, the magic lines of environment modules and zsh completion
# functions, and the mode line of a stylesheet, which sits outside its block
PRELUDE = re.compile(r"^(#!|<\?xml |(?i:<!doctype )|\s+\S|#%Module|#compdef |/\*\s+-\*- [a-z]+ -\*-\s+\*/$)")
# the mode line; the mode is in lower case, unless it carries variables, like the image name of
# a dockerfile, whose double marker it keeps
MODE = re.compile(r"^#?(?P<marker>.*?)-\*- (?P<mode>[a-z0-9+\-]+|[a-z\-]+: .+) -\*-$")

# a commit subject names the area it touches, then says what the commit does: an area is a
# path, a package, or a c++ namespace, and several can be listed separated by commas
AREA = r"[A-Za-z0-9_.{}/+\-]+(::[A-Za-z0-9_.{}/+\-]+)*"
# so a subject is one or more areas, a colon, and what the commit does
SUBJECT = re.compile(rf"^{AREA}(, ?{AREA})*: \S")
# the authors whose commits the commit check leaves alone: bots whose messages are beyond our
# control, such as the dependency updates dependabot proposes
BOTS = ("dependabot[bot]",)
# the copyright line, with any closing year
COPYRIGHT = re.compile(r"\(c\) 2010-(\d{4}) The PyLith authors and contributors$")


class ConventionError(Exception):
    """
    The base class of the reasons a check could not run
    """


class ToolError(ConventionError):
    """
    A tool a check relies on is missing
    """


def run(*, command: list, cwd: str) -> subprocess.CompletedProcess:
    """
    Run {command} in {cwd} and capture what it says, complaining if the tool is missing
    """
    # attempt to
    try:
        # run the command
        return subprocess.run(command, cwd=cwd, capture_output=True, text=True)
    # if the tool is not installed
    except FileNotFoundError as error:
        # say which one, and where to get it
        raise ToolError(f"'{command[0]}' is not installed; see etc/ci/conventions/requirements.txt") from error


def tracked(*, root: str, extensions: tuple = (), exempt: tuple = ()) -> list:
    """
    The files git tracks under {root}, narrowed to {extensions} and outside the {exempt} patterns
    """
    # ask git
    names = run(command=["git", "ls-files"], cwd=root).stdout.split("\n")
    # keep the ones of interest
    return [
        name
        for name in names
        if name
        and (not extensions or name.endswith(extensions))
        and not any(fnmatch.fnmatch(name, pattern) for pattern in exempt)
        and os.path.isfile(os.path.join(root, name))
    ]


def exemptions(*, root: str) -> tuple:
    """
    The patterns of the files the preamble and c++ formatting checks leave alone, as listed in
    the {exempt} table of {tool.conventions} in {pyproject.toml}
    """
    # the configuration of the repository
    path = os.path.join(root, "pyproject.toml")
    # a repository without one
    if not os.path.isfile(path):
        # exempts nothing
        return ()
    # read it
    with open(path, mode="rb") as stream:
        # all of it
        configuration = tomllib.load(stream)
    # hand off the patterns
    return tuple(configuration.get("tool", {}).get("conventions", {}).get("exempt", []))


def skips(*, root: str) -> tuple:
    """
    The patterns of the files the spell check leaves alone, as listed in {skip} of
    {tool.codespell} in {pyproject.toml}
    """
    # the configuration of the repository
    path = os.path.join(root, "pyproject.toml")
    # a repository without one
    if not os.path.isfile(path):
        # skips nothing
        return ()
    # read it
    with open(path, mode="rb") as stream:
        # all of it
        configuration = tomllib.load(stream)
    # the patterns are a single string, separated by commas
    patterns = configuration.get("tool", {}).get("codespell", {}).get("skip", "")
    # split them, and hand off the ones that are not blank
    return tuple(pattern.strip() for pattern in patterns.split(",") if pattern.strip())


def skipped(*, name: str, patterns: tuple) -> bool:
    """
    Check whether the tracked file {name} matches any of the skip {patterns}, the way codespell
    matches them when it walks a directory: against the name of the file and of each directory
    above it, and against each of their paths, with and without a leading './'
    """
    # the directories above the file, and the file itself
    parts = name.split("/")
    # the paths that lead to each of them
    paths = ["/".join(parts[: n + 1]) for n in range(len(parts))]
    # every spelling a pattern may be written against
    candidates = parts + paths + [f"./{path}" for path in paths]
    # the file is skipped when any pattern matches any of them
    return any(fnmatch.fnmatch(candidate, pattern) for candidate in candidates for pattern in patterns)


def spelling(*, root: str, base: str | None) -> list:
    """
    Spell check the prose of every tracked file, with the exclusions in {pyproject.toml}
    """
    # the files codespell would skip if it walked the tree; named on the command line, it
    # checks them regardless, so they are left out here
    patterns = skips(root=root)
    # the files
    files = [name for name in tracked(root=root) if not skipped(name=name, patterns=patterns)]
    # check them
    result = run(command=["codespell", "--toml", "pyproject.toml", *files], cwd=root)
    # every line of output is a misspelling
    problems = [line for line in result.stdout.split("\n") if line.strip()]
    # add the advice, if there is anything to advise about
    if problems:
        # a word codespell does not know is either a typo or belongs on the reviewed lists
        problems.append(
            "fix the typos; a real word goes in etc/ci/conventions/codespell-words.txt, and a "
            "line that is right as it is goes in etc/ci/conventions/codespell-lines.txt"
        )
    # hand off the complaints
    return problems


def python(*, root: str, base: str | None) -> list:
    """
    Check that every tracked python source and stub is formatted by black
    """
    # the files
    files = tracked(root=root, extensions=(".py", ".pyi"))
    # check them
    # the exclusions in {pyproject.toml} apply even to files named on the command line
    result = run(command=["black", *BLACK_OPTIONS, "--check", *files], cwd=root)
    # black names each file it would change
    problems = [line for line in result.stderr.split("\n") if line.startswith("would reformat")]
    # a failure that named no file is a failure of black itself
    if result.returncode != 0 and not problems:
        # so report what it said, or at least that it failed
        problems = [line for line in result.stderr.split("\n") if line.strip()] or [
            f"black exited with status {result.returncode}"
        ]
    # hand off the complaints
    return problems


def cxx(*, root: str, base: str | None) -> list:
    """
    Check that every tracked c++ source is formatted by clang-format
    """
    # the files
    files = tracked(root=root, extensions=CXX, exempt=exemptions(root=root))
    # check them
    result = run(command=["clang-format", "--dry-run", *files], cwd=root)
    # clang-format warns once per change it would make; one complaint per file is enough
    names = sorted({line.split(":")[0] for line in result.stderr.split("\n") if ": warning:" in line})
    # hand off the complaints
    return [f"{name}: not formatted by clang-format" for name in names]


def preambles(*, root: str, base: str | None) -> list:
    """
    Check that every tracked file opens with the standard preamble and closes with an end of
    file marker, and that the files touched since {base} carry the current year
    """
    # the complaints
    problems = []
    # the files, minus the exemptions of the repository
    files = tracked(root=root, exempt=exemptions(root=root))
    # the files touched since the base, when there is one
    touched = set(changed(root=root, base=base)) if base else set()
    # the current year
    year = str(datetime.date.today().year)
    # go through the files
    for name in files:
        # attempt to
        try:
            # read the file
            with open(os.path.join(root, name), encoding="utf-8") as stream:
                # all of it
                text = stream.read()
        # a file that is not text
        except UnicodeDecodeError:
            # has no preamble to check
            continue
        # an opaque format is checked only when the file carries a preamble anyway, as a
        # makefile named for the data it builds or a python source named for its role does
        if name.endswith(OPAQUE) and not any(COPYRIGHT_LINE.match(line) for line in text.split("\n")[:15]):
            # the rest are left alone
            continue
        # check its layout
        problems.extend(f"{name}: {problem}" for problem in layout(name=name, text=text))
        # the copyright line sits in the preamble, near the top
        copyright = next(
            (COPYRIGHT.search(line) for line in text.split("\n")[:15] if COPYRIGHT.search(line)),
            None,
        )
        # a file touched by the work
        if copyright is not None and name in touched and copyright.group(1) != year:
            # carries the current year
            problems.append(f"{name}: touched, but its copyright ends in {copyright.group(1)}")
    # hand off the complaints
    return problems


def layout(*, name: str, text: str) -> list:
    """
    Check the preamble and postamble of the file {name}, whose contents are {text}: the mode,
    coding, spacer, author and copyright lines, then the blank lines the kind of file calls for
    before the body and before the end of file marker
    """
    # an empty file
    if not text:
        # has no preamble to speak of
        return ["no preamble"]
    # the complaints
    problems = []
    # the file ends with a single newline
    if not text.endswith("\n") or text.endswith("\n\n"):
        # or else its last line is not where it should be
        problems.append("does not end with exactly one newline")
    # split into lines, without the empty one after the final newline
    lines = text.rstrip("\n").split("\n")
    # the copyright line sits in the preamble, near the top
    idx, copyright = next(
        ((i, match) for i, line in enumerate(lines[:15]) if (match := COPYRIGHT_LINE.match(line))),
        (None, None),
    )
    # a file without one
    if idx is None or copyright is None:
        # is missing its preamble, and nothing else can be said about it
        return problems + ["no 'Copyright (c) 2010-YYYY The PyLith authors and contributors' in its preamble"]
    # the comment marker is whatever precedes the copyright
    prefix = copyright.group("prefix")
    # a block comment preamble opens on a line of its own
    opener = next((i for i in range(idx - 1, -1, -1) if lines[i].strip() in OPENERS), None)
    # the comment lines start right after the opener, or with the first line that carries the
    # marker, past the lines that may precede the preamble
    start = opener + 1 if opener is not None else 0
    # skip the prelude of a line comment preamble
    while opener is None and start < idx and PRELUDE.match(lines[start]):
        # one line at a time
        start += 1
    # every line before the comment, other than the opener, belongs in the prelude
    prelude = [line for line in lines[: opener if opener is not None else start] if line]
    # so anything else
    if not all(PRELUDE.match(line) for line in prelude):
        # is out of place
        problems.append("unexpected lines above the preamble")
    # the comment lines, through the copyright
    comment = lines[start : idx + 1]
    # the marker, without its trailing space
    marker = prefix.rstrip()
    # the mode line can sit in the prelude of a stylesheet
    moded = any("-*-" in line for line in prelude)
    # otherwise, it opens the comment
    if not moded:
        # get it
        mode = MODE.match(comment[0]) if comment else None
        # it must be there, in lower case
        if not mode or not comment[0].startswith(marker):
            # or else the language of the file is not recorded
            problems.append("the preamble does not open with a lower case mode line")
        # move on
        comment = comment[1:]
    # the rest of the comment: coding, spacer, at least one line of text, a spacer, and the copyright
    expected = [f"{prefix}-*- coding: utf-8 -*-", marker]
    # compare
    if (
        comment[:2] != expected
        or len(comment) < 5
        or comment[-2] != marker
        or any(not line.strip() or line.rstrip() == marker for line in comment[2:-2])
    ):
        # the preamble has lines missing, extra, or out of order
        problems.append("the preamble is not mode, coding, spacer, text, spacer, copyright, and license")
    # past the copyright
    nxt = idx + 1
    # the license line follows the copyright
    if nxt < len(lines) and lines[nxt].startswith(prefix) and LICENSE_LINE.match(lines[nxt][len(prefix) :]):
        # so move past it
        nxt += 1
    # or else
    else:
        # the license is not recorded
        problems.append("no 'SPDX-License-Identifier: LICENSE' line after the copyright")
    # a block comment closes on the next line
    if opener is not None:
        # make sure it does
        if nxt >= len(lines) or lines[nxt].strip() not in CLOSERS:
            # or else the preamble runs on
            problems.append("the preamble block does not close after the license")
        # and move past it
        nxt += 1
    # count the blank lines after the preamble
    blanks = 0
    # by walking down
    while nxt + blanks < len(lines) and not lines[nxt + blanks]:
        # one at a time
        blanks += 1
    # the first line of the body
    body = lines[nxt + blanks] if nxt + blanks < len(lines) else ""
    # a c++ header keeps its code guard next to the preamble, as does a markdown document, and
    # black keeps the body of a python stub one line away from it
    gap = 1 if body.startswith("// code guard") or name.endswith((".md", ".pyi", *YAML)) else 2
    # check
    if blanks != gap:
        # and complain
        problems.append(f"{blanks} blank lines after the preamble, instead of {gap}")
    # the last line
    last = lines[-1]
    # is the closing marker
    if "End of file" not in last:
        # or else the file is not closed
        problems.append("does not end with an 'End of file' marker")
        # and its spacing means nothing
        return problems
    # count the blank lines before it
    blanks = 0
    # by walking up
    while blanks < len(lines) - 1 and not lines[-2 - blanks]:
        # one at a time
        blanks += 1
    # black puts a single blank line after imports, so python decides for itself, and a single
    # one at the end of a python stub
    allowed = (1, 2) if name.endswith(".py") else (1,) if name.endswith((".pyi", ".md", *YAML)) else (2,)
    # a file with no body shares the blank lines after the preamble
    if len(lines) - 1 - blanks == nxt:
        # so its count is the one that matters
        allowed = (blanks,)
    # check
    if blanks not in allowed:
        # and complain
        problems.append(f"{blanks} blank lines before the End of file marker, instead of {allowed[-1]}")
    # hand off the complaints
    return problems


def narrative(*, root: str, base: str | None) -> list:
    """
    Check that every statement of the python sources touched since {base} has a comment on the
    line above it, so the comments read as the specification of the code
    """
    # without a base there is nothing to check
    if not base:
        # so say nothing
        return []
    # the exemptions of the repository
    exempt = exemptions(root=root)
    # the complaints
    problems = []
    # go through the python sources touched by the work
    for name in changed(root=root, base=base):
        # skip the ones that are not python, or that the repository exempts
        if not name.endswith(".py") or any(fnmatch.fnmatch(name, pattern) for pattern in exempt):
            # by moving on
            continue
        # read the file
        with open(os.path.join(root, name), encoding="utf-8") as stream:
            # all of it
            text = stream.read()
        # attempt to
        try:
            # parse it
            tree = ast.parse(text, filename=name)
        # a file that does not parse
        except SyntaxError as error:
            # cannot be checked, which is a problem of its own
            problems.append(f"{name}:{error.lineno}: does not parse")
            # so move on
            continue
        # the lines that end with a comment of their own, as class attributes often do
        trailing = {
            token.start[0]
            for token in tokenize.generate_tokens(io.StringIO(text).readline)
            if token.type == tokenize.COMMENT and token.line.strip()[0] != "#"
        }
        # check its statements
        problems.extend(
            f"{name}:{line}: no comment above this statement"
            for line in uncommented(lines=text.split("\n"), trailing=trailing, body=tree.body)
        )
    # hand off the complaints
    return problems


def uncommented(*, lines: list, trailing: set, body: list) -> list:
    """
    The line numbers of the statements in {body}, and in the blocks nested in them, that have no
    comment on the line above them, or one of their own at the end of their first line; {lines}
    is the text of the file, and {trailing} the numbers of the lines that end with a comment
    """
    # nothing comes before the first statement of a block
    previous = None
    # the line numbers of the offenders
    offenders = []
    # go through the statements, remembering the one before
    for statement in body:
        # a docstring documents its block
        if isinstance(statement, ast.Expr) and isinstance(getattr(statement, "value", None), ast.Constant):
            # and needs no comment
            pass
        # a function or a class carries a docstring
        elif isinstance(statement, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            # instead of a comment
            pass
        # an import that follows another
        elif isinstance(statement, (ast.Import, ast.ImportFrom)) and isinstance(previous, (ast.Import, ast.ImportFrom)):
            # shares its comment, if the two are adjacent
            if not adjacent(statement=statement, previous=previous):
                # otherwise it starts a group that needs one of its own
                offenders.extend(missing(lines=lines, trailing=trailing, statement=statement))
        # the metadata of a declaration, such as the documentation of a trait, that follows it
        elif annotates(statement=statement, previous=previous):
            # shares its comment
            pass
        # everything else
        else:
            # needs a comment of its own
            offenders.extend(missing(lines=lines, trailing=trailing, statement=statement))
        # the blocks nested in the statement
        for field in ("body", "orelse", "finalbody"):
            # are checked the same way
            offenders.extend(uncommented(lines=lines, trailing=trailing, body=getattr(statement, field, None) or []))
        # and so are the handlers of a {try}
        for handler in getattr(statement, "handlers", None) or []:
            # each with its own block
            offenders.extend(uncommented(lines=lines, trailing=trailing, body=handler.body))
        # and the cases of a {match}
        for case in getattr(statement, "cases", None) or []:
            # each with its own block
            offenders.extend(uncommented(lines=lines, trailing=trailing, body=case.body))
        # remember this one
        previous = statement
    # hand off the offenders
    return offenders


def subject(*, statement: ast.stmt | None) -> str | None:
    """
    The name that {statement} declares, as in {t = pyre.properties.str()}, or annotates, as in
    {t.doc = "..."} or {t.aliases |= {"a"}}; {None} for any other statement
    """
    # an augmented assignment has a single target
    if isinstance(statement, ast.AugAssign):
        # which is the one to look at
        target = statement.target
    # so does a simple assignment that sets one name
    elif isinstance(statement, ast.Assign) and len(statement.targets) == 1:
        # get it
        target = statement.targets[0]
    # nothing else has a subject
    else:
        # so say so
        return None
    # a declaration names its subject
    if isinstance(target, ast.Name):
        # directly
        return target.id
    # an annotation sets an attribute of it
    if isinstance(target, ast.Attribute) and isinstance(target.value, ast.Name):
        # so the subject is the owner of the attribute
        return target.value.id
    # anything else has no subject
    return None


def annotates(*, statement: ast.stmt, previous: ast.stmt | None) -> bool:
    """
    Check whether {statement} sets an attribute of the subject of {previous}, on the line right
    after it, the way the metadata of a trait follows its declaration
    """
    # the statement must set an attribute, or add to one
    if isinstance(statement, ast.AugAssign):
        # in which case its single target
        target = statement.target
    # a simple assignment
    elif isinstance(statement, ast.Assign):
        # sets its first one
        target = statement.targets[0]
    # and nothing else
    else:
        # is an annotation
        return False
    # the target must be an attribute
    if not isinstance(target, ast.Attribute):
        # or it is not an annotation
        return False
    # of the subject of the statement right above it
    return (
        adjacent(statement=statement, previous=previous)
        and subject(statement=statement) is not None
        and subject(statement=statement) == subject(statement=previous)
    )


def adjacent(*, statement: ast.stmt, previous: ast.stmt | None) -> bool:
    """
    Check whether {statement} starts on the line right after the end of {previous}
    """
    # a statement with nothing before it, or after one whose extent is unknown
    if previous is None or previous.end_lineno is None:
        # has no neighbor
        return False
    # otherwise, the two are neighbors when no line separates them
    return statement.lineno - previous.end_lineno == 1


def missing(*, lines: list, trailing: set, statement: ast.stmt) -> list:
    """
    The line number of {statement}, if the line above it, or above its decorators, is not a
    comment, and its first line does not end with one either
    """
    # a statement that carries its comment at the end of its line
    if statement.lineno in trailing:
        # is commented
        return []
    # a decorated statement starts with its first decorator
    first = min([statement.lineno, *(d.lineno for d in getattr(statement, "decorator_list", []))])
    # the line above it, counting from one
    above = lines[first - 2].strip() if first > 1 else ""
    # is a comment, unless it is code that was commented out, or the statement is an offender
    return [] if above.startswith("#") and not disabled(comment=above) else [statement.lineno]


def disabled(*, comment: str) -> bool:
    """
    Check whether {comment} is a line of code that was commented out, rather than prose
    """
    # the text of the comment
    text = comment.lstrip("#").strip()
    # prose rarely carries the punctuation of a call, an assignment, or a subscript
    if not any(mark in text for mark in "(=["):
        # so this is prose
        return False
    # attempt to
    try:
        # parse it as python
        ast.parse(text)
    # if it does not parse
    except SyntaxError:
        # it is prose
        return False
    # otherwise, it is code
    return True


def commits(*, root: str, base: str | None) -> list:
    """
    Check the commits since {base}: a subject that names its area and says what the commit does,
    and nothing after it; the commits of the bots in {BOTS}, whose messages are beyond our
    control, are left out
    """
    # without a base there is nothing to check
    if not base:
        # so say nothing
        return []
    # the commits since the base, oldest first, as hash, author, subject, and body, merges left out
    log = run(
        command=[
            "git",
            "log",
            "--no-merges",
            "--reverse",
            "--format=%h%x00%an%x00%s%x00%b%x01",
            f"{base}..HEAD",
        ],
        cwd=root,
    ).stdout
    # the complaints
    problems = []
    # go through the commits
    for record in log.split("\x01"):
        # skip the separator after the last one
        if not record.strip():
            # by moving on
            continue
        # unpack
        sha, author, subject, body = record.strip("\n").split("\x00")
        # a bot writes its messages its own way, and that way may change without notice
        if author in BOTS:
            # so its commits are not held to the conventions
            continue
        # the subject names its area and says what the commit does
        if not SUBJECT.match(subject):
            # or else it does not read like the rest of the history
            problems.append(f"{sha}: '{subject}' is not of the form 'area: what the commit does'")
    # hand off the complaints
    return problems


def changed(*, root: str, base: str) -> list:
    """
    The files added or modified since {base}
    """
    # ask git
    result = run(command=["git", "diff", "--name-only", "--diff-filter=AM", f"{base}...HEAD"], cwd=root)
    # hand off the names
    return [name for name in result.stdout.split("\n") if name]


def parse() -> argparse.Namespace:
    """
    Read the command line
    """
    # make a parser
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    # the base of the work
    parser.add_argument(
        "--base",
        default=None,
        help="the branch the work is based on, e.g. origin/main; enables the commit and year checks",
    )
    # a subset of the checks
    parser.add_argument(
        "--only",
        default=",".join(CHECKS),
        help=f"a comma separated subset of the checks: {', '.join(CHECKS)}",
    )
    # parse and hand off
    return parser.parse_args()


def main() -> int:
    """
    Run the requested checks and report what each one found
    """
    # read the command line
    options = parse()
    # the top of the checkout
    root = run(command=["git", "rev-parse", "--show-toplevel"], cwd=os.getcwd()).stdout.strip()
    # the checks, by name
    table = {name: globals()[name] for name in CHECKS}
    # the requested ones
    requested = [name for name in options.only.split(",") if name]
    # a name that is not a check is a mistake
    unknown = [name for name in requested if name not in table]
    # so refuse it
    if unknown:
        # say which
        print(f"unknown checks: {', '.join(unknown)}; the checks are {', '.join(CHECKS)}")
        # and bail
        return 2
    # the number of checks that found something
    failed = 0
    # go through them
    for name in requested:
        # attempt to
        try:
            # run the check
            problems = table[name](root=root, base=options.base)
        # if it could not run
        except ConventionError as error:
            # it counts as a failure
            problems = [str(error)]
        # report the verdict
        print(f"{name}: {'ok' if not problems else f'{len(problems)} problems'}")
        # and the details
        for problem in problems:
            # one per line
            print(f"    {problem}")
        # tally
        failed += 1 if problems else 0
    # the exit status says whether everything passed
    return 1 if failed else 0


# entry point
if __name__ == "__main__":
    # run the checks and report their verdict
    sys.exit(main())


# End of file
