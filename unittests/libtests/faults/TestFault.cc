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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestFault.hh" // Implementation of class methods

#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

#include <string> // USES std::string

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFault );

// ----------------------------------------------------------------------
// Test id()
void
pylith::faults::TestFault::testID(void)
{ // testID
  const int id = 346;
  FaultCohesiveKin fault;
  fault.id(id);
  
  CPPUNIT_ASSERT_EQUAL(id, fault.id());
} // testID

// ----------------------------------------------------------------------
// Test label()
void
pylith::faults::TestFault::testLabel(void)
{ // testLabel
  const std::string label = "the_database";
  FaultCohesiveKin fault;
  fault.label(label.c_str());
  
  CPPUNIT_ASSERT_EQUAL(label, std::string(fault.label()));
} // testLabel
    

// End of file 
