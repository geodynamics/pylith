// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace utils {
        class TestEventLogger;
    }
}

class pylith::utils::TestEventLogger : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Test constructor.
    static
    void testConstructor(void);

    /// Test get/setClassName().
    static
    void testClassName(void);

    /// Test initialize().
    static
    void testInitialize(void);

    /// Test registerEvent().
    static
    void testRegisterEvent(void);

    /// Test getEventId().
    static
    void testGetEventId(void);

    /// Test eventBegin() and eventEnd().
    static
    void testEventLogging(void);

    /// Test registerStage().
    static
    void testRegisterStage(void);

    /// Test getStageId().
    static
    void testGetStageId(void);

    /// Test stagePush() and stagePop().
    static
    void testStageLogging(void);

}; // class TestEventLogging

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestEventLogger::testConstructor", "[TestEventLogger]") {
    pylith::utils::TestEventLogger::testConstructor();
}
TEST_CASE("TestEventLogger::testClassName", "[TestEventLogger]") {
    pylith::utils::TestEventLogger::testClassName();
}
TEST_CASE("TestEventLogger::testInitialize", "[TestEventLogger]") {
    pylith::utils::TestEventLogger::testInitialize();
}
TEST_CASE("TestEventLogger::testRegisterEvent", "[TestEventLogger]") {
    pylith::utils::TestEventLogger::testRegisterEvent();
}
TEST_CASE("TestEventLogger::testGetEventId", "[TestEventLogger]") {
    pylith::utils::TestEventLogger::testGetEventId();
}
TEST_CASE("TestEventLogger::testEventLogging", "[TestEventLogger]") {
    pylith::utils::TestEventLogger::testEventLogging();
}
TEST_CASE("TestEventLogger::testRegisterStage", "[TestEventLogger]") {
    pylith::utils::TestEventLogger::testRegisterStage();
}
TEST_CASE("TestEventLogger::testGetStageId", "[TestEventLogger]") {
    pylith::utils::TestEventLogger::testGetStageId();
}
TEST_CASE("TestEventLogger::testStageLogging", "[TestEventLogger]") {
    pylith::utils::TestEventLogger::testStageLogging();
}

// ------------------------------------------------------------------------------------------------
// Test constructor.
void
pylith::utils::TestEventLogger::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    EventLogger logger;

    PYLITH_METHOD_END;
} // testConstructor


// ------------------------------------------------------------------------------------------------
// Test get/setClassName().
void
pylith::utils::TestEventLogger::testClassName(void) {
    PYLITH_METHOD_BEGIN;

    EventLogger logger;
    CHECK(std::string("") == std::string(logger.getClassName()));

    const std::string& name = "my name";
    logger.setClassName(name.c_str());
    CHECK(name == std::string(logger.getClassName()));

    PYLITH_METHOD_END;
} // testClassName


// ------------------------------------------------------------------------------------------------
// Test initialize().
void
pylith::utils::TestEventLogger::testInitialize(void) {
    PYLITH_METHOD_BEGIN;

    EventLogger logger;
    logger.setClassName("my class");
    CHECK(0 == logger._classId);
    logger.initialize();
    CHECK(logger._classId);

    PYLITH_METHOD_END;
} // testInitialize


// ------------------------------------------------------------------------------------------------
// Test registerEvent().
void
pylith::utils::TestEventLogger::testRegisterEvent(void) {
    PYLITH_METHOD_BEGIN;

    EventLogger logger;
    logger.setClassName("my class");
    logger.initialize();

    const int numEvents = 3;
    const char* events[numEvents] = { "event A", "event B", "event C" };
    int ids[numEvents];

    for (int i = 0; i < numEvents; ++i) {
        ids[i] = logger.registerEvent(events[i]);
    }

    int i = 0;
    for (EventLogger::map_event_type::iterator e_iter = logger._events.begin(); e_iter != logger._events.end(); ++e_iter, ++i) {
        CHECK(std::string(events[i]) == e_iter->first);
        CHECK(ids[i] == e_iter->second);
    } // for

    PYLITH_METHOD_END;
} // testRegisterEvent


// ------------------------------------------------------------------------------------------------
// Test getEventId().
void
pylith::utils::TestEventLogger::testGetEventId(void) {
    PYLITH_METHOD_BEGIN;

    EventLogger logger;
    logger.setClassName("my class");
    logger.initialize();

    const int numEvents = 3;
    const char* events[numEvents] = { "event A", "event B", "event C" };

    for (int i = 0; i < numEvents; ++i) {
        logger.registerEvent(events[i]);
    }

    const int order[numEvents] = { 1, 0, 2 };
    int ids[numEvents];
    for (int i = 0; i < numEvents; ++i) {
        ids[order[i]] = logger.getEventId(events[order[i]]);
    }

    int i = 0;
    for (EventLogger::map_event_type::iterator e_iter = logger._events.begin(); e_iter != logger._events.end(); ++e_iter, ++i) {
        CHECK(e_iter->second == ids[i]);
    }

    PYLITH_METHOD_END;
} // testGetEventId


// ------------------------------------------------------------------------------------------------
// Test eventBegin() and eventEnd().
void
pylith::utils::TestEventLogger::testEventLogging(void) {
    PYLITH_METHOD_BEGIN;

    EventLogger logger;
    logger.setClassName("my class");
    logger.initialize();

    const int numEvents = 3;
    const char* events[numEvents] = { "event A", "event B", "event C" };
    int ids[numEvents];

    for (int i = 0; i < numEvents; ++i) {
        ids[i] = logger.registerEvent(events[i]);
    }

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


// ------------------------------------------------------------------------------------------------
// Test registerStage().
void
pylith::utils::TestEventLogger::testRegisterStage(void) {
    PYLITH_METHOD_BEGIN;

    EventLogger logger;
    logger.setClassName("my class");
    logger.initialize();

    const int numStages = 3;
    const char* stages[numStages] = { "stage A1", "stage B1", "stage C1" };
    int ids[numStages];

    for (int i = 0; i < numStages; ++i) {
        ids[i] = logger.registerStage(stages[i]);
    }

    int i = 0;
    for (EventLogger::map_event_type::iterator s_iter = logger._stages.begin(); s_iter != logger._stages.end(); ++s_iter, ++i) {
        CHECK(std::string(stages[i]) == s_iter->first);
        CHECK(ids[i] == s_iter->second);
    } // for

    PYLITH_METHOD_END;
} // testRegisterStage


// ------------------------------------------------------------------------------------------------
// Test stageId().
void
pylith::utils::TestEventLogger::testGetStageId(void) {
    PYLITH_METHOD_BEGIN;

    EventLogger logger;
    logger.setClassName("my class");
    logger.initialize();

    const int numStages = 3;
    const char* stages[numStages] = { "stage A2", "stage B2", "stage C2" };

    for (int i = 0; i < numStages; ++i) {
        logger.registerStage(stages[i]);
    }

    const int order[numStages] = { 1, 0, 2 };
    int ids[numStages];
    for (int i = 0; i < numStages; ++i) {
        ids[order[i]] = logger.getStageId(stages[order[i]]);
    }

    int i = 0;
    for (EventLogger::map_event_type::iterator s_iter = logger._stages.begin(); s_iter != logger._stages.end(); ++s_iter, ++i) {
        CHECK(s_iter->second == ids[i]);
    }

    const int idNew = logger.getStageId("stage D2");
    CHECK(idNew == logger.getStageId("stage D2"));

    PYLITH_METHOD_END;
} // testStageId


// ------------------------------------------------------------------------------------------------
// Test statePush() and statePop().
void
pylith::utils::TestEventLogger::testStageLogging(void) {
    PYLITH_METHOD_BEGIN;

    EventLogger logger;
    logger.setClassName("my class");
    logger.initialize();

    const int numStages = 3;
    const char* stages[numStages] = { "stage A3", "stage B3", "stage C3" };
    int ids[numStages];

    for (int i = 0; i < numStages; ++i) {
        ids[i] = logger.registerStage(stages[i]);
    }

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
