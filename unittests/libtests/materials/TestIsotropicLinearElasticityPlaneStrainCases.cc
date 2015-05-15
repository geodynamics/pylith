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

#include "TestIsotropicLinearElasticityPlaneStrain.hh" // Implementation of cases

#include "data/IsotropicLinearElasticityPlaneStrainData_Tri3.hh"

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {

    class TestIsotropicLinearElasticityPlaneStrain_Tri3 : public TestIsotropicLinearElasticityPlaneStrain
    { // TestIsotropicLinearElasticityPlaneStrain_Tri3

      CPPUNIT_TEST_SUB_SUITE( TestIsotropicLinearElasticityPlaneStrain_Tri3,  TestIsotropicLinearElasticityPlaneStrain);
      CPPUNIT_TEST_SUITE_END();

      void setUp(void) {
	TestIsotropicLinearElasticityPlaneStrain::setUp();
	_data = new IsotropicLinearElasticityPlaneStrainData_Tri3();
	_initializeMin();
      } // setUp

    }; // TestIsotropicLinearElasticityPlaneStrain_Tri3
    CPPUNIT_TEST_SUITE_REGISTRATION( TestIsotropicLinearElasticityPlaneStrain_Tri3 );


  } // materials
} // pylith


// End of file 
