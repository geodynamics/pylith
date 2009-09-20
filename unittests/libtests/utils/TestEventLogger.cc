// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestEventLogger.hh" // Implementation of class methods

#include "pylith/utils/EventLogger.hh" // USES EventLogger


// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::utils::TestEventLogger );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::utils::TestEventLogger::testConstructor(void)
{ // testConstructor
  EventLogger logger;
} // testConstructor

// ----------------------------------------------------------------------
// Test className().
void
pylith::utils::TestEventLogger::testClassName(void)
{ // testClassName
  EventLogger logger;
  CPPUNIT_ASSERT_EQUAL(std::string(""), std::string(logger.className()));

  const std::string& name = "my name";
  logger.className(name.c_str());
  CPPUNIT_ASSERT_EQUAL(name, std::string(logger.className()));
} // testClassName

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::utils::TestEventLogger::testInitialize(void)
{ // testInitialize
  EventLogger logger;
  logger.className("my class");
  CPPUNIT_ASSERT_EQUAL(0, logger._classId);
  logger.initialize();
  CPPUNIT_ASSERT(0 != logger._classId);
} // testInitialize

// ----------------------------------------------------------------------
// Test registerEvent().
void
pylith::utils::TestEventLogger::testRegisterEvent(void)
{ // testRegisterEvent
  EventLogger logger;
  logger.className("my class");
  logger.initialize();

  const char* events[] = { "event A", "event B", "event C" };
  const int numEvents = 3;
  int ids[numEvents];

  for (int i=0; i < numEvents; ++i)
    ids[i] = logger.registerEvent(events[i]);

  int i = 0;
  for (EventLogger::map_event_type::iterator e_iter=logger._events.begin();
       e_iter != logger._events.end();
       ++e_iter, ++i) {
    CPPUNIT_ASSERT_EQUAL(std::string(events[i]), e_iter->first);
    CPPUNIT_ASSERT_EQUAL(ids[i], e_iter->second);
  } // for
} // testRegisterEvent

// ----------------------------------------------------------------------
// Test eventId().
void
pylith::utils::TestEventLogger::testEventId(void)
{ // testEventId
  EventLogger logger;
  logger.className("my class");
  logger.initialize();

  const char* events[] = { "event A", "event B", "event C" };
  const int numEvents = 3;

  for (int i=0; i < numEvents; ++i)
    logger.registerEvent(events[i]);

  const int order[] = { 1, 0, 2 };
  int ids[numEvents];
  for (int i=0; i < numEvents; ++i)
    ids[order[i]] = logger.eventId(events[order[i]]);

  int i = 0;
  for (EventLogger::map_event_type::iterator e_iter=logger._events.begin();
       e_iter != logger._events.end();
       ++e_iter, ++i)
    CPPUNIT_ASSERT_EQUAL(e_iter->second, ids[i]);
} // testEventId

// ----------------------------------------------------------------------
// Test eventBegin() and eventEnd().
void
pylith::utils::TestEventLogger::testEventLogging(void)
{ // testEventLogging
  EventLogger logger;
  logger.className("my class");
  logger.initialize();

  const char* events[] = { "event A", "event B", "event C" };
  const int numEvents = 3;
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
} // testEventLogging

// ----------------------------------------------------------------------
// Test registerStage().
void
pylith::utils::TestEventLogger::testRegisterStage(void)
{ // testRegisterStage
  EventLogger logger;
  logger.className("my class");
  logger.initialize();

  const char* stages[] = { "stage A1", "stage B1", "stage C1" };
  const int numStages = 3;
  int ids[numStages];

  for (int i=0; i < numStages; ++i)
    ids[i] = logger.registerStage(stages[i]);

  int i = 0;
  for (EventLogger::map_event_type::iterator s_iter=logger._stages.begin();
       s_iter != logger._stages.end();
       ++s_iter, ++i) {
    CPPUNIT_ASSERT_EQUAL(std::string(stages[i]), s_iter->first);
    CPPUNIT_ASSERT_EQUAL(ids[i], s_iter->second);
  } // for
} // testRegisterStage

// ----------------------------------------------------------------------
// Test stageId().
void
pylith::utils::TestEventLogger::testStageId(void)
{ // testStageId
  EventLogger logger;
  logger.className("my class");
  logger.initialize();

  const char* stages[] = { "stage A2", "stage B2", "stage C2" };
  const int numStages = 3;

  for (int i=0; i < numStages; ++i)
    logger.registerStage(stages[i]);


  const int order[] = { 1, 0, 2 };
  int ids[numStages];
  for (int i=0; i < numStages; ++i)
    ids[order[i]] = logger.stageId(stages[order[i]]);

  int i = 0;
  for (EventLogger::map_event_type::iterator s_iter=logger._stages.begin();
       s_iter != logger._stages.end();
       ++s_iter, ++i)
    CPPUNIT_ASSERT_EQUAL(s_iter->second, ids[i]);
  
  const int idNew = logger.stageId("stage D2");
  CPPUNIT_ASSERT_EQUAL(idNew, logger.stageId("stage D2"));
} // testStageId

// ----------------------------------------------------------------------
// Test statePush() and statePop().
void
pylith::utils::TestEventLogger::testStageLogging(void)
{ // testStageLogging
  EventLogger logger;
  logger.className("my class");
  logger.initialize();

  const char* stages[] = { "stage A3", "stage B3", "stage C3" };
  const int numStages = 3;
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
} // testStageLogging


// End of file 
