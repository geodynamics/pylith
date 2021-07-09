#!/usr/bin/env python
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ======================================================================
#

# @file tests/pytests/utils/TestEventLogger.py

# @brief Unit testing of EventLogger object.

import unittest


# ----------------------------------------------------------------------
class TestEventLogger(unittest.TestCase):
    """Unit testing of EventLogger object.
    """

    def test_constructor(self):
        from pylith.utils.EventLogger import EventLogger
        logger = EventLogger()
        return

    def test_className(self):
        from pylith.utils.EventLogger import EventLogger
        logger = EventLogger()
        name = "my class"
        logger.setClassName(name)
        self.assertEqual(name, logger.getClassName())
        return

    def test_initialize(self):
        from pylith.utils.EventLogger import EventLogger
        logger = EventLogger()
        logger.setClassName("logging A")
        logger.initialize()
        return

    def test_registerEvent(self):
        from pylith.utils.EventLogger import EventLogger
        logger = EventLogger()
        logger.setClassName("logging A")
        logger.initialize()

        events = ["event 1", "event 2", "event 3"]
        id = {}
        for event in events:
            id[event] = logger.registerEvent(event)
        for event in events:
            self.assertEqual(id[event], logger.getEventId(event))
        return

    def test_eventLogging(self):
        from pylith.utils.EventLogger import EventLogger
        logger = EventLogger()
        logger.setClassName("logging A")
        logger.initialize()
        events = ["event 1", "event 2", "event 3"]
        for event in events:
            logger.registerEvent(event)

        logger.eventBegin("event 2")
        logger.eventEnd("event 2")

        logger.eventBegin("event 1")
        logger.eventEnd("event 1")

        logger.eventBegin("event 3")
        logger.eventBegin("event 1")
        logger.eventEnd("event 1")
        logger.eventEnd("event 3")
        return

    def test_registerStage(self):
        from pylith.utils.EventLogger import EventLogger
        logger = EventLogger()
        logger.setClassName("logging A")
        logger.initialize()

        stages = ["stage 1a", "stage 2a", "stage 3a"]
        id = {}
        for stage in stages:
            id[stage] = logger.registerStage(stage)
        for stage in stages:
            self.assertEqual(id[stage], logger.getStageId(stage))
        return

    def test_stageLogging(self):
        from pylith.utils.EventLogger import EventLogger
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
        return


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestEventLogger))

    from pylith.utils.PetscManager import PetscManager
    petsc = PetscManager()
    petsc.initialize()

    success = unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful()

    petsc.finalize()


# End of file
