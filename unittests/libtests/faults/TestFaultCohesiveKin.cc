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

#include "TestFaultCohesiveKin.hh" // Implementation of class methods

#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKin );

// ----------------------------------------------------------------------
// Test clone()
void
pylith::faults::TestFaultCohesiveKin::testClone(void)
{ // testClone
  const int id = 65;
  const std::string label = "fault ABC";

  FaultCohesiveKin fault;
  fault.id(id);
  fault.label(label.c_str());
  
  Fault* faultCopy = fault.clone();
  CPPUNIT_ASSERT(0 != faultCopy);
  
  CPPUNIT_ASSERT_EQUAL(id, faultCopy->id());
  CPPUNIT_ASSERT_EQUAL(label, faultCopy->label());
} // testClone


// End of file 
