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

#include "pylith/faults/EqKinSrc.hh" // USES EqKinSrc

#include <stdexcept> // TEMPORARY
// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKin );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestFaultCohesiveKin::testConstructor(void)
{ // testConstructor
  FaultCohesiveKin fault;
} // testConstructor

// ----------------------------------------------------------------------
// Test eqsrc().
void
pylith::faults::TestFaultCohesiveKin::testEqsrc(void)
{ // testEqsrc
  FaultCohesiveKin fault;

  EqKinSrc eqsrc;
  fault.eqsrc(&eqsrc);
  CPPUNIT_ASSERT(&eqsrc == fault._eqsrc);
} // testEqsrc

// ----------------------------------------------------------------------
// Test useLagrangeConstraints().
void
pylith::faults::TestFaultCohesiveKin::testUseLagrangeConstraints(void)
{ // testUseLagrangeConstraints
  FaultCohesiveKin fault;
  CPPUNIT_ASSERT_EQUAL(true, fault._useLagrangeConstraints());
} // testUseLagrangeConstraints

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestFaultCohesiveKin::testInitialize(void)
{ // testInitialize
  throw std::logic_error("Unit test not implemented.");
} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::faults::TestFaultCohesiveKin::testIntegrateResidual(void)
{ // testIntegrateResidual
  throw std::logic_error("Unit test not implemented.");
} // testIntegrateResidual

// ----------------------------------------------------------------------
// Test integrateJacobian().
void
pylith::faults::TestFaultCohesiveKin::testIntegrateJacobian(void)
{ // testIntegrateJacobian
  throw std::logic_error("Unit test not implemented.");
} // testIntegrateJacobian

// ----------------------------------------------------------------------
// Test setConstraintSizes().
void
pylith::faults::TestFaultCohesiveKin::testSetConstraintSizes(void)
{ // testSetConstraintSizes
  // Make sure the fiber dimension at each point is equal to the
  // spatial dimension of the mesh, because the Lagrange multiplier
  // formation does not eliminate any DOF from the system of
  // equations.

  throw std::logic_error("Unit test not implemented.");
} // testSetConstraintSizes

// ----------------------------------------------------------------------
// Test setField().
void
pylith::faults::TestFaultCohesiveKin::testSetField(void)
{ // testSetField
  throw std::logic_error("Unit test not implemented.");
} // testSetField

// ----------------------------------------------------------------------
// Initialize FaultCohesiveKin interface condition.
void
pylith::faults::TestFaultCohesiveKin::_initialize(ALE::Obj<ALE::Mesh>* mesh,
					FaultCohesiveKin* const fault) const
{ // _initialize
} // _initialize


// End of file 
