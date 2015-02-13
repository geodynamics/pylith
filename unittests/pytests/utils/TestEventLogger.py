#!/usr/bin/env python
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/pytests/utils/TestEventLogger.py

## @brief Unit testing of EventLogger object.

import unittest


# ----------------------------------------------------------------------
class TestEventLogger(unittest.TestCase):
  """
  Unit testing of EventLogger object.
  """
  

  def test_constructor(self):
    """
    Test constructor.
    """
    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    return


  def test_className(self):
    """
    Test className().
    """
    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    name = "my class"
    logger.className(name)
    self.assertEqual(name, logger.className())
    return


  def test_initialize(self):
    """
    Test initialize().
    """
    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("logging A")
    logger.initialize()
    return


  def test_registerEvent(self):
    """
    Test registerEvent().
    """
    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("logging A")
    logger.initialize()

    events = ["event 1" , "event 2" , "event 3"]
    id = {}
    for event in events:
      id[event] = logger.registerEvent(event)
    for event in events:
      self.assertEqual(id[event], logger.eventId(event))
    return


  def test_eventLogging(self):
    """
    Test eventBegin() and eventEnd().
    """
    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("logging A")
    logger.initialize()
    events = ["event 1" , "event 2" , "event 3"]
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
    """
    Test registerStage().
    """
    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("logging A")
    logger.initialize()

    stages = ["stage 1a" , "stage 2a" , "stage 3a"]
    id = {}
    for stage in stages:
      id[stage] = logger.registerStage(stage)
    for stage in stages:
      self.assertEqual(id[stage], logger.stageId(stage))
    return


  def test_stageLogging(self):
    """
    Test stagePush() and stagePop().
    """
    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("logging A")
    logger.initialize()
    stages = ["stage 1b" , "stage 2b" , "stage 3b"]
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


# End of file 
