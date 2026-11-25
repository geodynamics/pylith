// -*- c++ -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


#include <portinfo>

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/petsc/EventLogger.hh" // USES EventLogger

#include "pylith/exceptions/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"
#include <stdexcept> // USES std::runtime_error

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace petsc {
        class TestEventLogger;
    }
} // namespace pylith

class pylith::petsc::TestEventLogger : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
  public:
    /// Test constructor.
    static void testConstructor(void);

    /// Test get/setClassName().
    static void testClassName(void);

    /// Test initialize().
    static void testInitialize(void);

    /// Test registerEvent().
    static void testRegisterEvent(void);

    /// Test EventLogger::Scope.
    static void testEventScope(void);

    /// Test registerStage().
    static void testRegisterStage(void);

    /// Test getStageId().
    static void testGetStageId(void);

    /// Test stagePush() and stagePop().
    static void testStageLogging(void);

}; // class TestEventLogging

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestEventLogger::testConstructor",
          "[TestEventLogger]") {
    pylith::petsc::TestEventLogger::testConstructor();
}
TEST_CASE("TestEventLogger::testClassName",
          "[TestEventLogger]") {
    pylith::petsc::TestEventLogger::testClassName();
}
TEST_CASE("TestEventLogger::testInitialize",
          "[TestEventLogger]") {
    pylith::petsc::TestEventLogger::testInitialize();
}
TEST_CASE("TestEventLogger::testRegisterEvent",
          "[TestEventLogger]") {
    pylith::petsc::TestEventLogger::testRegisterEvent();
}
TEST_CASE("TestEventLogger::testEventScope",
          "[TestEventLogger]") {
    pylith::petsc::TestEventLogger::testEventScope();
}
TEST_CASE("TestEventLogger::testRegisterStage",
          "[TestEventLogger]") {
    pylith::petsc::TestEventLogger::testRegisterStage();
}
TEST_CASE("TestEventLogger::testGetStageId",
          "[TestEventLogger]") {
    pylith::petsc::TestEventLogger::testGetStageId();
}
TEST_CASE("TestEventLogger::testStageLogging",
          "[TestEventLogger]") {
    pylith::petsc::TestEventLogger::testStageLogging();
}

// ------------------------------------------------------------------------------------------------
// Test constructor.
void
pylith::petsc::TestEventLogger::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    EventLogger logger;

    PYLITH_METHOD_END;
} // testConstructor


// ------------------------------------------------------------------------------------------------
// Test get/setClassName().
void
pylith::petsc::TestEventLogger::testClassName(void) {
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
pylith::petsc::TestEventLogger::testInitialize(void) {
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
pylith::petsc::TestEventLogger::testRegisterEvent(void) {
    PYLITH_METHOD_BEGIN;

    EventLogger logger;
    logger.setClassName("my class");
    logger.initialize();

    const int numEvents = 3;
    const char* events[numEvents] = {"event A", "event B", "event C"};
    PetscLogEvent ids[numEvents];

    for (int i = 0; i < numEvents; ++i) {
        ids[i] = logger.registerEvent(events[i]);
    }

    int i = 0;
    for (EventLogger::map_event_type::iterator e_iter = logger._events.begin(); e_iter != logger._events.end();
         ++e_iter, ++i) {
        CHECK(std::string(events[i]) == e_iter->first);
        CHECK(ids[i] == e_iter->second);
    } // for

    PYLITH_METHOD_END;
} // testRegisterEvent


// ------------------------------------------------------------------------------------------------
// Test EventLogger::Scope.
void
pylith::petsc::TestEventLogger::testEventScope(void) {
    PYLITH_METHOD_BEGIN;

    EventLogger logger;
    logger.setClassName("my class");
    logger.initialize();
    const PetscLogEvent eventA = logger.registerEvent("scope event A");
    const PetscLogEvent eventB = logger.registerEvent("scope event B");

    { // Nested scopes.
        EventLogger::Scope scopeA(eventA);
        EventLogger::Scope scopeB(eventB);
    }

    try { // Scope ends event when an exception is thrown.
        EventLogger::Scope scopeA(eventA);
        throw std::runtime_error("test");
    } catch (const std::runtime_error&) {
    }

    // Scope ends event on early return.
    const auto earlyReturn = [&](void) {
        EventLogger::Scope scopeB(eventB);
        return;
    };
    earlyReturn();

    PYLITH_METHOD_END;
} // testEventScope


// ------------------------------------------------------------------------------------------------
// Test registerStage().
void
pylith::petsc::TestEventLogger::testRegisterStage(void) {
    PYLITH_METHOD_BEGIN;

    EventLogger logger;
    logger.setClassName("my class");
    logger.initialize();

    const int numStages = 3;
    const char* stages[numStages] = {"stage A1", "stage B1", "stage C1"};
    PetscLogStage ids[numStages];

    for (int i = 0; i < numStages; ++i) {
        ids[i] = logger.registerStage(stages[i]);
    }

    int i = 0;
    for (EventLogger::map_stage_type::iterator s_iter = logger._stages.begin(); s_iter != logger._stages.end();
         ++s_iter, ++i) {
        CHECK(std::string(stages[i]) == s_iter->first);
        CHECK(ids[i] == s_iter->second);
    } // for

    PYLITH_METHOD_END;
} // testRegisterStage


// ------------------------------------------------------------------------------------------------
// Test stageId().
void
pylith::petsc::TestEventLogger::testGetStageId(void) {
    PYLITH_METHOD_BEGIN;

    EventLogger logger;
    logger.setClassName("my class");
    logger.initialize();

    const int numStages = 3;
    const char* stages[numStages] = {"stage A2", "stage B2", "stage C2"};

    for (int i = 0; i < numStages; ++i) {
        logger.registerStage(stages[i]);
    }

    const int order[numStages] = {1, 0, 2};
    PetscLogStage ids[numStages];
    for (int i = 0; i < numStages; ++i) {
        ids[order[i]] = logger.getStageId(stages[order[i]]);
    }

    int i = 0;
    for (EventLogger::map_stage_type::iterator s_iter = logger._stages.begin(); s_iter != logger._stages.end();
         ++s_iter, ++i) {
        CHECK(s_iter->second == ids[i]);
    }

    const int idNew = logger.getStageId("stage D2");
    CHECK(idNew == logger.getStageId("stage D2"));

    PYLITH_METHOD_END;
} // testStageId


// ------------------------------------------------------------------------------------------------
// Test statePush() and statePop().
void
pylith::petsc::TestEventLogger::testStageLogging(void) {
    PYLITH_METHOD_BEGIN;

    EventLogger logger;
    logger.setClassName("my class");
    logger.initialize();

    const int numStages = 3;
    const char* stages[numStages] = {"stage A3", "stage B3", "stage C3"};
    PetscLogStage ids[numStages];

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
