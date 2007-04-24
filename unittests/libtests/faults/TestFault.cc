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

#include "TestFault.hh" // Implementation of class methods

#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

#include <string> // USES std::string
#include <stdexcept> // USES std::logic_error TEMPORARY

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
  
  CPPUNIT_ASSERT(id == fault.id());
} // testID

// ----------------------------------------------------------------------
// Test label()
void
pylith::faults::TestFault::testLabel(void)
{ // testLabel
  const std::string label = "the_database";
  FaultCohesiveKin fault;
  fault.label(label.c_str());
  
  CPPUNIT_ASSERT_EQUAL(label, fault.label());
} // testLabel
    
#if 0
// ----------------------------------------------------------------------
// Test initialize()
void
pylith::faults::TestFault::testInitialize(void)
{ // testInitialize
  throw std::logic_error("Need to implement test.");
} // testInitialize
#endif


// End of file 
