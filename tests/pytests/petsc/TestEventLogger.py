# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

import unittest

from pylith.testing.TestCases import make_suite
from pylith.petsc.EventLogger import EventLogger


class TestEventLogger(unittest.TestCase):
    """Unit testing of EventLogger object."""

    def test_constructor(self):
        logger = EventLogger()

    def test_className(self):
        logger = EventLogger()
        name = "my class"
        logger.setClassName(name)
        self.assertEqual(name, logger.getClassName())

    def test_initialize(self):
        logger = EventLogger()
        logger.setClassName("logging A")
        logger.initialize()

    def test_registerStage(self):
        logger = EventLogger()
        logger.setClassName("logging A")
        logger.initialize()

        stages = ["stage 1a", "stage 2a", "stage 3a"]
        id = {}
        for stage in stages:
            id[stage] = logger.registerStage(stage)
        for stage in stages:
            self.assertEqual(id[stage], logger.getStageId(stage))

    def test_stageLogging(self):
        logger = EventLogger()
        logger.setClassName("logging A")
        logger.initialize()
        stages = ["stage 1b", "stage 2b", "stage 3b"]
        for stage in stages:
            logger.registerStage(stage)

        logger.stagePush("stage 2b")
        logger.stagePop()

        logger.stagePush("stage 1b")
        logger.stagePop()

        logger.stagePush("stage 3b")
        logger.stagePush("stage 1b")
        logger.stagePop()
        logger.stagePop()


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [TestEventLogger]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    from pylith.petsc.Manager import Manager

    petsc = Manager()
    petsc.initialize()

    unittest.main(verbosity=2)

    petsc.finalize()


# End of file
