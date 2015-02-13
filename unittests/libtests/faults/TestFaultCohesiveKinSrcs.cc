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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestFaultCohesiveKinSrcs.hh" // Implementation of class methods

#include "pylith/faults/EqKinSrc.hh" // USES EqKinSrc
#include "pylith/faults/BruneSlipFn.hh" // USES BruneSlipFn

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveKinSrcs );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveKinSrcs::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestFaultCohesiveKin::setUp();
  _data = 0;

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

  PYLITH_METHOD_END;
} // setUp


// End of file 
