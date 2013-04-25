// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestEventLogger.hh" // Implementation of class methods

#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include "pylith/utils/petscerror.h" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::utils::TestEventLogger );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::utils::TestEventLogger::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  EventLogger logger;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test className().
void
pylith::utils::TestEventLogger::testClassName(void)
{ // testClassName
  PYLITH_METHOD_BEGIN;

  EventLogger logger;
  CPPUNIT_ASSERT_EQUAL(std::string(""), std::string(logger.className()));

  const std::string& name = "my name";
  logger.className(name.c_str());
  CPPUNIT_ASSERT_EQUAL(name, std::string(logger.className()));

  PYLITH_METHOD_END;
} // testClassName

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::utils::TestEventLogger::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;

  EventLogger logger;
  logger.className("my class");
  CPPUNIT_ASSERT_EQUAL(0, logger._classId);
  logger.initialize();
  CPPUNIT_ASSERT(logger._classId);

  PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
// Test registerEvent().
void
pylith::utils::TestEventLogger::testRegisterEvent(void)
{ // testRegisterEvent
  PYLITH_METHOD_BEGIN;

  EventLogger logger;
  logger.className("my class");
  logger.initialize();

  const int numEvents = 3;
  const char* events[numEvents] = { "event A", "event B", "event C" };
  int ids[numEvents];

  for (int i=0; i < numEvents; ++i)
    ids[i] = logger.registerEvent(events[i]);

  int i = 0;
  for (EventLogger::map_event_type::iterator e_iter=logger._events.begin(); e_iter != logger._events.end(); ++e_iter, ++i) {
    CPPUNIT_ASSERT_EQUAL(std::string(events[i]), e_iter->first);
    CPPUNIT_ASSERT_EQUAL(ids[i], e_iter->second);
  } // for

  PYLITH_METHOD_END;
} // testRegisterEvent

// ----------------------------------------------------------------------
// Test eventId().
void
pylith::utils::TestEventLogger::testEventId(void)
{ // testEventId
  PYLITH_METHOD_BEGIN;

  EventLogger logger;
  logger.className("my class");
  logger.initialize();

  const int numEvents = 3;
  const char* events[numEvents] = { "event A", "event B", "event C" };

  for (int i=0; i < numEvents; ++i)
    logger.registerEvent(events[i]);

  const int order[numEvents] = { 1, 0, 2 };
  int ids[numEvents];
  for (int i=0; i < numEvents; ++i)
    ids[order[i]] = logger.eventId(events[order[i]]);

  int i = 0;
  for (EventLogger::map_event_type::iterator e_iter=logger._events.begin(); e_iter != logger._events.end(); ++e_iter, ++i)
    CPPUNIT_ASSERT_EQUAL(e_iter->second, ids[i]);

  PYLITH_METHOD_END;
} // testEventId

// ----------------------------------------------------------------------
// Test eventBegin() and eventEnd().
void
pylith::utils::TestEventLogger::testEventLogging(void)
{ // testEventLogging
  PYLITH_METHOD_BEGIN;

  EventLogger logger;
  logger.className("my class");
  logger.initialize();

  const int numEvents = 3;
  const char* events[numEvents] = { "event A", "event B", "event C" };
  int ids[numEvents];

  for (int i=0; i < numEvents; ++i)
    ids[i] = logger.registerEvent(events[i]);

  int event = ids[1];
  logger.eventBegin(event);
  logger.eventEnd(event);

  event = ids[0];
  logger.eventBegin(event);
  logger.eventEnd(event);

  event = ids[2];
  logger.eventBegin(event);
  int event2 = ids[0];
  logger.eventBegin(event2);
  logger.eventEnd(event2);
  logger.eventEnd(event);

  PYLITH_METHOD_END;
} // testEventLogging

// ----------------------------------------------------------------------
// Test registerStage().
void
pylith::utils::TestEventLogger::testRegisterStage(void)
{ // testRegisterStage
  PYLITH_METHOD_BEGIN;

  EventLogger logger;
  logger.className("my class");
  logger.initialize();

  const int numStages = 3;
  const char* stages[numStages] = { "stage A1", "stage B1", "stage C1" };
  int ids[numStages];

  for (int i=0; i < numStages; ++i)
    ids[i] = logger.registerStage(stages[i]);

  int i = 0;
  for (EventLogger::map_event_type::iterator s_iter=logger._stages.begin(); s_iter != logger._stages.end(); ++s_iter, ++i) {
    CPPUNIT_ASSERT_EQUAL(std::string(stages[i]), s_iter->first);
    CPPUNIT_ASSERT_EQUAL(ids[i], s_iter->second);
  } // for

  PYLITH_METHOD_END;
} // testRegisterStage

// ----------------------------------------------------------------------
// Test stageId().
void
pylith::utils::TestEventLogger::testStageId(void)
{ // testStageId
  PYLITH_METHOD_BEGIN;

  EventLogger logger;
  logger.className("my class");
  logger.initialize();

  const int numStages = 3;
  const char* stages[numStages] = { "stage A2", "stage B2", "stage C2" };

  for (int i=0; i < numStages; ++i)
    logger.registerStage(stages[i]);


  const int order[numStages] = { 1, 0, 2 };
  int ids[numStages];
  for (int i=0; i < numStages; ++i)
    ids[order[i]] = logger.stageId(stages[order[i]]);

  int i = 0;
  for (EventLogger::map_event_type::iterator s_iter=logger._stages.begin(); s_iter != logger._stages.end();  ++s_iter, ++i)
    CPPUNIT_ASSERT_EQUAL(s_iter->second, ids[i]);
  
  const int idNew = logger.stageId("stage D2");
  CPPUNIT_ASSERT_EQUAL(idNew, logger.stageId("stage D2"));

  PYLITH_METHOD_END;
} // testStageId

// ----------------------------------------------------------------------
// Test statePush() and statePop().
void
pylith::utils::TestEventLogger::testStageLogging(void)
{ // testStageLogging
  PYLITH_METHOD_BEGIN;

  EventLogger logger;
  logger.className("my class");
  logger.initialize();

  const int numStages = 3;
  const char* stages[numStages] = { "stage A3", "stage B3", "stage C3" };
  int ids[numStages];

  for (int i=0; i < numStages; ++i)
    ids[i] = logger.registerStage(stages[i]);

  int stage = ids[1];
  logger.stagePush(stage);
  logger.stagePop();

  stage = ids[0];
  logger.stagePush(stage);
  logger.stagePop();

  stage = ids[2];
  logger.stagePush(stage);
  int stage2 = ids[0];
  logger.stagePush(stage2);
  logger.stagePop();
  logger.stagePop();

  PYLITH_METHOD_END;
} // testStageLogging


// End of file 
