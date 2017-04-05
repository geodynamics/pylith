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

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::DiscretizeInfo

// ----------------------------------------------------------------------
namespace pylith {
    namespace materials {

        class TestIsotropicLinearElasticityPlaneStrain_UniformStrain : public TestIsotropicLinearElasticityPlaneStrain
        { // TestIsotropicLinearElasticityPlaneStrain_Tri3
protected:

            void setUp(void) {
                TestIsotropicLinearElasticityPlaneStrain::setUp();
                _mydata = new TestIsotropicLinearElasticityPlaneStrain_Data(); CPPUNIT_ASSERT(_mydata);

                _mydata->numSolnFields = 1;
                _mydata->useInertia = false;
                _mydata->useBodyForce = false;
                _mydata->useGravity = false;
                _mydata->useReferenceState = false;

                _mydata->materialLabel = "IsotropicLinearElascitity";
                _mydata->materialId = 24;
                _mydata->boundaryLabel = "boundary";

                _mydata->auxDBFilename = "data/IsotropicLinearElasticityPlaneStrain_UniStrain_aux.spatialdb";
                _mydata->solnDBFilename = "data/IsotropicLinearElasticityPlaneStrain_UniStrain_soln.spatialdb";
                _mydata->pertDBFilename = "data/IsotropicLinearElasticityPlaneStrain_UniStrain_pert.spatialdb";

                _mydata->lengthScale = 1.0e+03;
                _mydata->timeScale = 2.0;
                _mydata->densityScale = 3.0e+3;
                _mydata->pressureScale = 2.25e+10;

                _mydata->t = 1.0;
                _mydata->dt = 0.05;
                _mydata->tshift = 1.0 / _mydata->dt;

                _mydata->numAuxFields = 3;
                static const char* _auxFields[3] = {"density", "shear_modulus", "bulk_modulus"};
                _mydata->auxFields = _auxFields;
                static const pylith::topology::Field::DiscretizeInfo _auxDiscretizations[3] = {
                    {0, 1, true}, // density
                    {0, 1, true}, // shear_modulus
                    {0, 1, true}, // bulk_modulus
                };
                _mydata->auxDiscretizations = const_cast<pylith::topology::Field::DiscretizeInfo*>(_auxDiscretizations);

            } // setUp

        }; // TestIsotropicLinearElasticityPlaneStrain_UniformStrain


    } // materials
} // pylith


// ----------------------------------------------------------------------
namespace pylith {
    namespace materials {

        class TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP1 : public TestIsotropicLinearElasticityPlaneStrain_UniformStrain
        { // TestIsotropicLinearElasticityPlaneStrain_UniformStrain_Tri3

            CPPUNIT_TEST_SUB_SUITE( TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP1,  TestIsotropicLinearElasticityPlaneStrain);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestIsotropicLinearElasticityPlaneStrain_UniformStrain::setUp();
                CPPUNIT_ASSERT(_mydata);

                _mydata->meshFilename = "data/tri3_small.mesh";

                static const pylith::topology::Field::DiscretizeInfo _solnDiscretizations[1] = {
                    {1, 1, true}, // disp
                };
                _mydata->solnDiscretizations = const_cast<pylith::topology::Field::DiscretizeInfo*>(_solnDiscretizations);

                _initializeMin();
            } // setUp

        }; // TestIsotropicLinearElasticityPlaneStrain_Tri3
        CPPUNIT_TEST_SUITE_REGISTRATION( TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP1 );


    } // materials
} // pylith


// ----------------------------------------------------------------------
namespace pylith {
    namespace materials {

        class TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP2 : public TestIsotropicLinearElasticityPlaneStrain_UniformStrain
        { // TestIsotropicLinearElasticityPlaneStrain_UniformStrain_Tri3

            CPPUNIT_TEST_SUB_SUITE( TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP2,  TestIsotropicLinearElasticityPlaneStrain);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestIsotropicLinearElasticityPlaneStrain_UniformStrain::setUp();
                CPPUNIT_ASSERT(_mydata);

                _mydata->meshFilename = "data/tri3_small.mesh";
                static const pylith::topology::Field::DiscretizeInfo _solnDiscretizations[1] = {
                    {2, 2, true}, // disp
                };
                _mydata->solnDiscretizations = const_cast<pylith::topology::Field::DiscretizeInfo*>(_solnDiscretizations);

                static const pylith::topology::Field::DiscretizeInfo _auxDiscretizations[3] = {
                    {0, 2, true}, // density
                    {0, 2, true}, // shear_modulus
                    {0, 2, true}, // bulk_modulus
                };
                _mydata->auxDiscretizations = const_cast<pylith::topology::Field::DiscretizeInfo*>(_auxDiscretizations);

                _initializeMin();
            } // setUp

        }; // TestIsotropicLinearElasticityPlaneStrain_Tri3
        CPPUNIT_TEST_SUITE_REGISTRATION( TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP2 );


    } // materials
} // pylith


// End of file
