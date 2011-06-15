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

#include "TestBoundaryCondition.hh" // Implementation of class methods

#include "pylith/bc/DirichletBC.hh" // USES DirichletBC

#include <string> // USES std::string

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestBoundaryCondition );

// ----------------------------------------------------------------------
// Test label().
void
pylith::bc::TestBoundaryCondition::testLabel(void)
{ // testLabel
  const std::string& label = "the_database";
  DirichletBC bc;
  bc.label(label.c_str());
  
  CPPUNIT_ASSERT_EQUAL(label, std::string(bc.label()));
} // testLabel
    

// End of file 
