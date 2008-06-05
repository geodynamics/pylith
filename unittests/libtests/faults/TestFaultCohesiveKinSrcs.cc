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

#include "TestFaultCohesiveKinSrcs.hh" // Implementation of class methods

#include "pylith/faults/EqKinSrc.hh" // USES EqKinSrc
#include "pylith/faults/BruneSlipFn.hh" // USES BruneSlipFn

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinSrcs );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinSrcs::setUp(void)
{ // setUp
  TestFaultCohesiveKin::setUp();
  _data = 0;
  _quadrature = 0;

  // Cleanup previous allocation
  int nsrcs = _eqsrcs.size();
  for (int i=0; i < nsrcs; ++i)
    delete _eqsrcs[i];
  nsrcs = _slipfns.size();
  for (int i=0; i < nsrcs; ++i)
    delete _slipfns[i];

  nsrcs = 2;
  _eqsrcs.resize(nsrcs);
  _eqsrcs[0] = new EqKinSrc();
  _eqsrcs[0]->originTime(0.5);
  _eqsrcs[1] = new EqKinSrc();
  _eqsrcs[1]->originTime(0.0);

  _slipfns.resize(nsrcs);
  _slipfns[0] = new BruneSlipFn();
  _slipfns[1] = new BruneSlipFn();
} // setUp


// End of file 
