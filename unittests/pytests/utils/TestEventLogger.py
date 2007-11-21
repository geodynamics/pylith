#!/usr/bin/env python
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
#
# ======================================================================
#

## @file unittests/pytests/utils/TestPescManager.py

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
    Test setClassName() and getClassName.
    """
    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    name = "my class"
    logger.setClassName(name)
    self.assertEqual(name, logger.getClassName())
    return


  def test_initialize(self):
    """
    Test initialize().
    """
    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.setClassName("logging A")
    logger.initialize()
    return


  def test_registerEvent(self):
    """
    Test registerEvent().
    """
    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.setClassName("logging A")
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
    logger.setClassName("logging A")
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


# End of file 
