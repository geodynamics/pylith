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
import datetime
import os
import shutil
import subprocess
import sys
import tempfile

# the checker under test, next to this file
CHECKER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "check.py")
# the c++ formatting rules of the repository, which the scratch repositories borrow
CLANG_FORMAT = os.path.join(os.path.dirname(CHECKER), "..", "..", "..", ".clang-format")
# the current year, which the preamble of a touched file must carry
YEAR = datetime.date.today().year
# a misspelling for the spell checker to catch, assembled here so the source that plants it
# does not trip the spell check of this repository
MISSPELLING = "t" + "eh"


class SelfTestError(Exception):
    """
    The checker did not report what it should have
    """


def preamble(*, comment: str, language: str) -> str:
    """
    The standard preamble of a source file, in the given {comment} syntax
    """
    # assemble it line by line
    return (
        f"{comment} -*- {language} -*-\n"
        f"{comment} -*- coding: utf-8 -*-\n"
        f"{comment}\n"
        f"{comment} This code is part of PyLith.\n"
        f"{comment}\n"
        f"{comment} Copyright (c) 2010-{YEAR} The PyLith authors and contributors\n"
        f"{comment} SPDX-License-Identifier: MIT\n"
    )


def git(*args: str, cwd: str) -> str:
    """
    Run git in {cwd} and hand back what it says
    """
    # run it, failing loudly, since a scratch repository that cannot be built tests nothing
    result = subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True, check=True)
    # hand off its output
    return result.stdout


def write(*, root: str, name: str, text: str) -> None:
    """
    Write {text} to the file {name} under {root}
    """
    # make room for it
    os.makedirs(os.path.dirname(os.path.join(root, name)) or root, exist_ok=True)
    # and write it
    with open(os.path.join(root, name), mode="w", encoding="utf-8") as stream:
        # all at once
        stream.write(text)
    # all done
    return


def scaffold(*, root: str) -> None:
    """
    Lay out a scratch repository whose single commit satisfies every convention
    """
    # make the repository
    git("init", "--quiet", "--initial-branch=main", cwd=root)
    # with an identity for its commits: a name
    git("config", "user.name", "self test", cwd=root)
    # and an address
    git("config", "user.email", "selftest@example.com", cwd=root)
    # the configuration the checks read
    write(
        root=root,
        name="pyproject.toml",
        text=preamble(comment="#", language="toml")
        + '\n\n[tool.black]\nline-length = 100\n\n[tool.codespell]\nignore-words = "words.txt"\n'
        + '\n[tool.conventions]\nexempt = ["vendor/*"]\n'
        + "\n\n# End of file\n",
    )
    # an empty list of known good words
    write(root=root, name="words.txt", text="")
    # the c++ formatting rules
    shutil.copy(CLANG_FORMAT, os.path.join(root, ".clang-format"))
    # a python source that follows the conventions
    write(
        root=root,
        name="good.py",
        text=preamble(comment="#", language="python")
        + '\n\ndef answer():\n    """\n    The answer\n    """\n    # hand it off\n    return 42\n'
        + "\n\n# End of file\n",
    )
    # and a c++ one
    write(
        root=root,
        name="good.cc",
        text=preamble(comment="//", language="c++")
        + "\n\n// the answer\nint\nanswer() {\n    // hand it off\n    return 42;\n}\n\n\n// End of file\n",
    )
    # stage them
    git("add", "pyproject.toml", "words.txt", ".clang-format", "good.py", "good.cc", cwd=root)
    # and commit them
    git("commit", "--quiet", "-m", "base: the starting point", cwd=root)
    # all done
    return


def check(*, root: str, base: str | None = None) -> tuple:
    """
    Run the checker in {root} and hand back its exit status and report
    """
    # the command line
    command = [sys.executable, CHECKER] + (["--base", base] if base else [])
    # run it
    result = subprocess.run(command, cwd=root, capture_output=True, text=True)
    # hand off the verdict
    return result.returncode, result.stdout


def expect(*, report: str, check: str, needle: str) -> None:
    """
    Make sure the {check} section of {report} mentions {needle}
    """
    # find the section of the check
    lines = report.split("\n")
    # its header
    start = next((i for i, line in enumerate(lines) if line.startswith(f"{check}:")), None)
    # a check that did not report at all
    if start is None:
        # is a failure
        raise SelfTestError(f"'{check}' did not report:\n{report}")
    # its details are the indented lines that follow
    details = []
    # collect them
    for line in lines[start + 1 :]:
        # until the next section
        if not line.startswith("    "):
            # which ends this one
            break
        # otherwise, keep it
        details.append(line)
    # the needle has to be among them
    if not any(needle in line for line in details):
        # or else the check missed what it was supposed to catch
        raise SelfTestError(f"'{check}' did not mention '{needle}':\n{report}")
    # all done
    return


def clean() -> None:
    """
    A repository that follows the conventions passes every check
    """
    # in a scratch directory
    with tempfile.TemporaryDirectory() as root:
        # lay out the repository
        scaffold(root=root)
        # and check it
        status, report = check(root=root, base="main")
        # everything passed
        if status != 0:
            # or else the checker complains about good work
            raise SelfTestError(f"a clean repository failed:\n{report}")
    # all done
    return


def violations() -> None:
    """
    Each kind of violation is caught by the check that looks for it
    """
    # in a scratch directory
    with tempfile.TemporaryDirectory() as root:
        # lay out the repository
        scaffold(root=root)
        # work on a branch
        git("switch", "--quiet", "-c", "work", cwd=root)
        # a python source with a misspelling, bad formatting, and no closing marker
        write(
            root=root,
            name="bad.py",
            text=preamble(comment="#", language="python") + f"\n\n# {MISSPELLING} answer\nx=[1,2 ,3]\n",
        )
        # a c++ source that is not formatted
        write(
            root=root,
            name="bad.cc",
            text=preamble(comment="//", language="c++") + "\nint answer(){return 42;}\n\n\n// End of file\n",
        )
        # an existing file touched without bringing its year up to date
        with open(os.path.join(root, "good.cc"), encoding="utf-8") as stream:
            # read it
            text = stream.read()
        # age its copyright, and change something in it
        write(
            root=root,
            name="good.cc",
            text=text.replace(f"2010-{YEAR}", "2010-2000").replace("the answer", "the reply"),
        )
        # stage it all
        git("add", "bad.py", "bad.cc", "good.cc", cwd=root)
        # and commit it with a subject that names no area, and a body
        git("commit", "--quiet", "-m", "fixed some things", "-m", "a body", cwd=root)
        # check
        status, report = check(root=root, base="main")
        # the checker failed
        if status == 0:
            # or else it let everything through
            raise SelfTestError(f"violations passed:\n{report}")
        # and every check caught its violation: the misspelling
        expect(report=report, check="spelling", needle=MISSPELLING)
        # the python formatting
        expect(report=report, check="python", needle="bad.py")
        # the c++ formatting
        expect(report=report, check="cxx", needle="bad.cc")
        # the missing closing marker
        expect(report=report, check="preambles", needle="bad.py: does not end")
        # the stale copyright year
        expect(report=report, check="preambles", needle="good.cc: touched")
        # the subject that names no area
        expect(report=report, check="commits", needle="is not of the form")
    # all done
    return


def layouts() -> None:
    """
    Each departure from the standard preamble and postamble layout is caught, and the files the
    repository exempts are left alone
    """
    # the departures, by file, with what the checker says about each
    cases = {
        # a mode in mixed case
        "upper.py": (
            preamble(comment="#", language="Python") + "\n\nx = 1\n\n\n# End of file\n",
            "lower case mode line",
        ),
        # no coding line
        "coding.sh": (
            "# -*- bash -*-\n#\n# text\n#\n# Copyright (c) 2010-2026 The PyLith authors and contributors\n# SPDX-License-Identifier: MIT\n"
            + "\n\necho\n\n\n# End of file\n",
            "mode, coding, spacer, text, spacer, copyright, and license",
        ),
        # two blank lines after the preamble of a yaml file, whose formatter keeps one
        "loose.yaml": (
            preamble(comment="#", language="yaml") + "\n\nkey: value\n\n\n# End of file\n",
            "2 blank lines after the preamble, instead of 1",
        ),
        # a single blank line before the closing marker of a c++ source
        "short.cc": (
            preamble(comment="//", language="c++") + "\n\nint x = 0;\n\n// End of file\n",
            "1 blank lines before the End of file marker",
        ),
        # two blank lines between the preamble and the code guard of a header
        "guard.h": (
            preamble(comment="//", language="c++") + "\n\n// code guard\n#pragma once\n\n\n// End of file\n",
            "2 blank lines after the preamble, instead of 1",
        ),
        # a makefile named for the data it builds, whose extension is otherwise left alone
        "build.data": (
            preamble(comment="#", language="makefile") + "\nall:\n\n\n# End of file\n",
            "1 blank lines after the preamble",
        ),
        # a python stub spaced like a python source
        "loose.pyi": (
            preamble(comment="#", language="python") + "\n\nx: int\n\n\n# End of file\n",
            "2 blank lines after the preamble, instead of 1",
        ),
        # extra blank lines at the end
        "trail.toml": (
            preamble(comment="#", language="toml") + "\n\nkey = 1\n\n\n# End of file\n\n",
            "exactly one newline",
        ),
    }
    # the departures the checker must not flag: a header with its code guard, a markdown
    # document, a python stub as black lays it out, a dockerfile, a script, a stylesheet, and a
    # file the repository exempts
    fine = {
        "tight.yaml": preamble(comment="#", language="yaml") + "\nkey: value\n\n# End of file\n",
        "tight.pyi": preamble(comment="#", language="python") + "\nx: int\n\n# End of file\n",
        "guard.icc": preamble(comment="//", language="c++") + "\n// code guard\n#pragma once\n\n\n// End of file\n",
        "doc.md": "<!--\n-*- markdown -*-\n-*- coding: utf-8 -*-\n\n"
        + "text\n\n"
        + f"Copyright (c) 2010-{YEAR} The PyLith authors and contributors\nSPDX-License-Identifier: MIT\n-->\n\n# title\n\n<!-- End of file -->\n",
        "Dockerfile": '## -*- docker-image-name: "x:y" -*-\n'
        + preamble(comment="#", language="x").split("\n", 1)[1]
        + "\n\nFROM ubuntu\n\n\n# End of file\n",
        "run.sh": "#!/bin/bash\n" + preamble(comment="#", language="bash") + "\n\necho\n\n\n# End of file\n",
        "look.css": "/*  -*- css -*-  */\n/*\n * -*- coding: utf-8 -*-\n *\n"
        + " * text\n *\n"
        + f" * Copyright (c) 2010-{YEAR} The PyLith authors and contributors\n * SPDX-License-Identifier: MIT\n */\n\n\na {{\n}}\n\n\n/* End of file */\n",
        "_complete": "#compdef complete\n"
        + preamble(comment="#", language="zsh")
        + "\n\n_arguments\n\n\n# End of file\n",
        "page.html": "<!DOCTYPE html>\n\n<!--\n-*- html -*-\n-*- coding: utf-8 -*-\n\n"
        + "text\n\n"
        + f"Copyright (c) 2010-{YEAR} The PyLith authors and contributors\nSPDX-License-Identifier: MIT\n-->\n\n\n<html>\n</html>\n\n\n<!-- End of file -->\n",
        "vendor/theirs.js": "// no preamble here\n",
        "values.csv": "1,2,3\n",
    }
    # in a scratch directory
    with tempfile.TemporaryDirectory() as root:
        # lay out the repository
        scaffold(root=root)
        # add the files
        for name, text in [(name, text) for name, (text, _) in cases.items()] + list(fine.items()):
            # one by one
            write(root=root, name=name, text=text)
        # and let git know about them
        git("add", *cases, *fine, cwd=root)
        # check the preambles
        result = subprocess.run(
            [sys.executable, CHECKER, "--only", "preambles"],
            cwd=root,
            capture_output=True,
            text=True,
        )
        # every departure was caught
        for name, (_, needle) in cases.items():
            # with the right complaint on its own line
            if not any(f"{name}: " in line and needle in line for line in result.stdout.split("\n")):
                # or else the checker missed it
                raise SelfTestError(f"'{name}' was not flagged for '{needle}':\n{result.stdout}")
        # and nothing else was
        for name in fine:
            # by name
            if f"{name}: " in result.stdout:
                # or else the checker complains about good work
                raise SelfTestError(f"'{name}' was flagged:\n{result.stdout}")
    # all done
    return


def narratives() -> None:
    """
    A statement without a comment above it is caught in the python sources the work touches,
    and the forms the conventions allow are left alone
    """
    # in a scratch directory
    with tempfile.TemporaryDirectory() as root:
        # lay out the repository
        scaffold(root=root)
        # work on a branch
        git("switch", "--quiet", "-c", "work", cwd=root)
        # the body of a python source, one line per entry, with what the check must say
        body = [
            ("# support", None),
            ("import os", None),
            ("import sys", None),
            ("", None),
            ("", None),
            ("class gadget:", None),
            ('    """a gadget"""', None),
            ("", None),
            ("    # the size of the gadget", None),
            ("    size = os.sep", None),
            ('    size.doc = "how big"', None),
            ("    color = sys.platform  # the color of the gadget", None),
            ("", None),
            ("", None),
            ("def run():", None),
            ('    """', None),
            ("    Run", None),
            ('    """', None),
            ("    # print(os.sep)", None),
            ("    value = 1", "flag"),
            ("    total = value + 1", "flag"),
            ("    # all done", None),
            ("    return total", None),
        ]
        # the lines of the body start right after the preamble and the two blank lines
        first = len(preamble(comment="#", language="python").split("\n")) + 2
        # write the source
        write(
            root=root,
            name="story.py",
            text=preamble(comment="#", language="python")
            + "\n\n"
            + "\n".join(line for line, _ in body)
            + "\n\n\n# End of file\n",
        )
        # stage it
        git("add", "story.py", cwd=root)
        # and commit it
        git("commit", "--quiet", "-m", "story: a source with gaps in its narrative", cwd=root)
        # check
        status, report = check(root=root, base="main")
        # go through the body
        for offset, (line, verdict) in enumerate(body):
            # the statements the check must catch
            if verdict == "flag":
                # are named by their line
                expect(report=report, check="narrative", needle=f"story.py:{first + offset}:")
            # and the rest
            elif f"story.py:{first + offset}:" in report:
                # are not
                raise SelfTestError(f"'narrative' flagged '{line.strip()}':\n{report}")
    # all done
    return


def bots() -> None:
    """
    The commits of a bot are not held to the commit conventions, and the same message from a
    person still is
    """
    # in a scratch directory
    with tempfile.TemporaryDirectory() as root:
        # lay out the repository
        scaffold(root=root)
        # work on a branch
        git("switch", "--quiet", "-c", "work", cwd=root)
        # the message a dependency update comes with: no area, and a body
        message = ["-m", "Bump black from 26.5.1 to 26.6.0", "-m", "the release notes"]
        # a bot makes a change
        write(root=root, name="words.txt", text="tread\n")
        # stages it
        git("add", "words.txt", cwd=root)
        # and commits it with that message
        git("-c", "user.name=dependabot[bot]", "commit", "--quiet", *message, cwd=root)
        # which passes
        status, report = check(root=root, base="main")
        # every check
        if status != 0:
            # or else the bot is held to conventions it cannot follow
            raise SelfTestError(f"a bot commit failed:\n{report}")
        # a person makes another change
        write(root=root, name="words.txt", text="tread\nraison\n")
        # stages it
        git("add", "words.txt", cwd=root)
        # and commits it with the same message
        git("commit", "--quiet", *message, cwd=root)
        # which does not pass
        status, report = check(root=root, base="main")
        # the commit check catches it
        expect(report=report, check="commits", needle="is not of the form")
    # all done
    return


def skips() -> None:
    """
    The spell check leaves alone the files its skip patterns match, by the name of the file or of
    a directory above it, or by their paths, and still checks the rest
    """
    # the files with misspellings that the patterns cover, by the pattern that covers each
    covered = {
        # the name of the file, anywhere in the tree
        "web/package-lock.json": "package-lock.json",
        # the name of a directory above it
        "src/generated/notes.txt": "generated",
        # its path, the way the codespell documentation writes it
        "share/data/notes.txt": "./share/data/*",
    }
    # in a scratch directory
    with tempfile.TemporaryDirectory() as root:
        # lay out the repository
        scaffold(root=root)
        # skip the files the patterns cover
        with open(os.path.join(root, "pyproject.toml"), encoding="utf-8") as stream:
            # read the configuration
            text = stream.read()
        # and add the patterns to it
        write(
            root=root,
            name="pyproject.toml",
            text=text.replace(
                'ignore-words = "words.txt"\n',
                f'ignore-words = "words.txt"\nskip = "{",".join(covered.values())}"\n',
            ),
        )
        # plant a misspelling in each of the covered files
        for name in covered:
            # one by one
            write(root=root, name=name, text=f"{MISSPELLING}\n")
        # and in one that no pattern covers
        write(root=root, name="src/notes.txt", text=f"{MISSPELLING}\n")
        # let git know about them
        git("add", "pyproject.toml", *covered, "src/notes.txt", cwd=root)
        # check the spelling
        result = subprocess.run(
            [sys.executable, CHECKER, "--only", "spelling"],
            cwd=root,
            capture_output=True,
            text=True,
        )
        # the file no pattern covers was checked
        expect(report=result.stdout, check="spelling", needle="src/notes.txt")
        # and the covered ones were not
        for name in covered:
            # by name
            if f"{name}:" in result.stdout:
                # or else the skip patterns were ignored
                raise SelfTestError(f"'{name}' was not skipped:\n{result.stdout}")
        # the files the checker reported, from the details of its spelling section
        checked = {line.strip().split(":")[0] for line in result.stdout.split("\n") if "==>" in line}
        # every file in the repository is tracked, so codespell walking it applies the same
        # patterns itself, and must report the same files; a release of codespell that matches
        # its patterns differently fails here, rather than in the spell check of a repository
        walk = subprocess.run(
            ["codespell", "--toml", "pyproject.toml", "."],
            cwd=root,
            capture_output=True,
            text=True,
        )
        # the files the walk reported, without the leading './' of its paths
        walked = {line.split(":")[0].removeprefix("./") for line in walk.stdout.split("\n") if "==>" in line}
        # the two must agree
        if checked != walked:
            # or else the checker no longer skips what codespell skips
            raise SelfTestError(f"the checker reported {sorted(checked)}, codespell walking the tree {sorted(walked)}")
    # all done
    return


def main() -> int:
    """
    Run the self tests
    """
    # attempt to
    try:
        # check a clean repository
        clean()
        # one full of violations
        violations()
        # one with every kind of preamble departure
        layouts()
        # one with commits by a bot
        bots()
        # and one with files the spell check skips
        skips()
    # if the checker misbehaved
    except SelfTestError as error:
        # say how
        print(f"selftest: {error}")
        # and fail
        return 1
    # all good
    return 0


# entry point
if __name__ == "__main__":
    # run the self tests and report their verdict
    sys.exit(main())


# End of file
