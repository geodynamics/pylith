# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

import unittest

from pylith.testing.TestCases import make_suite
from pylith.utils.utils import PylithVersion

class TestPylithVersion(unittest.TestCase):

    def test_isRelease(self):
        isRelease = PylithVersion.isRelease()


    def test_version(self):
        version = PylithVersion.version()
        # Check that version is of the form X.X.X
        import re
        match = re.search("[0-9]+\\.[0-9]+\\.[0-9]+", version)
        self.assertFalse(match is None)


    def test_gitVersion(self):
        revision = PylithVersion.gitRevision()
        if PylithVersion.isRelease():
            self.assertEqual("unknown", revision)
        else:
            # Check that revision is of the form v2.1.3-16-g9323114
            import re
            match = re.search("v[0-9]+\\.[0-9]+\\.[0-9]+", revision)
            if match is None:
                match = re.search("v[0-9]+\\.[0-9]+\\.[0-9]+-[0-9]+-g[0-9,a-z]+", revision)
            self.assertFalse(match is None)


    def test_gitHash(self):
        tag = PylithVersion.gitHash()
        if PylithVersion.isRelease():
            self.assertEqual("unknown", tag)
        else:
            # Check form of hash
            import re
            match = re.search("[0-9,a-z]+", tag)
            self.assertFalse(match is None)


    def test_gitDate(self):
        value = PylithVersion.gitDate()
        if PylithVersion.isRelease():
            self.assertEqual("unknown", value)
        else:
            # Check form of datetime
            import datetime
            fields = value.split()
            d = datetime.datetime.strptime(fields[0], "%Y-%m-%d")
            t = datetime.datetime.strptime(fields[1], "%H:%M:%S")


    def test_gitBranch(self):
        branch = PylithVersion.gitBranch()
        if PylithVersion.isRelease():
            self.assertEqual("unknown", branch)
        else:
            self.assertFalse(len(branch) == 0)


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestPylithVersion]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file 
