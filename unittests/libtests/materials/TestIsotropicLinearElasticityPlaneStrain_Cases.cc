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

#include "TestIsotropicLinearElasticityPlaneStrain_Data.hh"

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::DiscretizeInfo

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {

    class TestIsotropicLinearElasticityPlaneStrain_Tri3 : public TestIsotropicLinearElasticityPlaneStrain
    { // TestIsotropicLinearElasticityPlaneStrain_Tri3

      CPPUNIT_TEST_SUB_SUITE( TestIsotropicLinearElasticityPlaneStrain_Tri3,  TestIsotropicLinearElasticityPlaneStrain);
      CPPUNIT_TEST_SUITE_END();

      void setUp(void) {
	const bool useInertia = false;
	const bool useBodyForce = false;
	const bool useInitialState = false;
	
	TestIsotropicLinearElasticityPlaneStrain::setUp();
	_data = new TestIsotropicLinearElasticityPlaneStrain_Data(useInertia, useBodyForce, useInitialState);


	_data->meshFilename = "data/tri3_small.mesh";
	_data-> materialLabel = "IsotropicLinearElascitity";
	_data->materialId = 24;
	_data->boundaryLabel = "boundary";

	_data->auxDBFilename = "data/IsotropicLinearElasticityPlaneStrain_UniStrain_aux.spatialdb";
	_data->solnDBFilename = "data/IsotropicLinearElasticityPlaneStrain_UniStrain_soln.spatialdb";
	_data->pertDBFilename = "data/IsotropicLinearElasticityPlaneStrain_UniStrain_pert.spatialdb";

	_data->lengthScale = 1.0e+03;
	_data->timeScale = 2.0;
	_data->densityScale = 3.0e+3;
	_data->pressureScale = 2.25e+10;

	_data->t = 1.0;
	_data->dt = 0.05;
	_data->tshift = 1.0/0.05;

	static const pylith::topology::Field::DiscretizeInfo _solnDiscretizations[2] = {
	  {1, 1, true}, // disp
	  {1, 1, true}, // vel
	};
	_data->solnDiscretizations = const_cast<pylith::topology::Field::DiscretizeInfo*>(_solnDiscretizations);

	_data->numAuxFields = 3;
	static const char* _auxFields[3] = {"density", "shear_modulus", "bulk_modulus"};
	static const pylith::topology::Field::DiscretizeInfo _auxDiscretizations[3] = {
	  {0, 1, true}, // density
	  {0, 1, true}, // shear_modulus
	  {0, 1, true}, // bulk_modulus
	};
	_data->auxFields = _auxFields;
	_data->auxDiscretizations = const_cast<pylith::topology::Field::DiscretizeInfo*>(_auxDiscretizations);

	_initializeMin();
      } // setUp

    }; // TestIsotropicLinearElasticityPlaneStrain_Tri3
    CPPUNIT_TEST_SUITE_REGISTRATION( TestIsotropicLinearElasticityPlaneStrain_Tri3 );


  } // materials
} // pylith


// End of file 
