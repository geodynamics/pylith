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

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
namespace pylith {
    namespace materials {

        class TestIsotropicLinearElasticityPlaneStrain_UniformStrain : public TestIsotropicLinearElasticityPlaneStrain {
protected:

            void setUp(void) {
                TestIsotropicLinearElasticityPlaneStrain::setUp();
                _mydata = new TestIsotropicLinearElasticityPlaneStrain_Data(); CPPUNIT_ASSERT(_mydata);

                // dimension set in base class.
                _mydata->useInertia = false;
                _mydata->useBodyForce = false;
                _mydata->useGravity = false;
                _mydata->useReferenceState = false;

                // meshFilename set in derived class.
                _mydata->materialLabel = "IsotropicLinearElascitity";
                _mydata->materialId = 24;
                _mydata->boundaryLabel = "boundary";

                CPPUNIT_ASSERT(_mydata->normalizer);
                _mydata->normalizer->lengthScale(1.0e+03);
                _mydata->normalizer->timeScale(2.0);
                _mydata->normalizer->densityScale(3.0e+3);
                _mydata->normalizer->pressureScale(2.25e+10);

                _mydata->t = 1.0;
                _mydata->dt = 0.05;
                _mydata->tshift = 1.0 / _mydata->dt;

                // solnDiscretizations set in derived class.
                _mydata->solnDBFilename = "data/IsotropicLinearElasticityPlaneStrain_UniformStrain_soln.spatialdb";
                _mydata->pertDBFilename = "data/IsotropicLinearElasticityPlaneStrain_UniformStrain_pert.spatialdb";

                _mydata->numAuxFields = 3;
                static const char* _auxFields[3] = {"density", "shear_modulus", "bulk_modulus"};
                _mydata->auxFields = _auxFields;
                static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // density
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // shear_modulus
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // bulk_modulus
                };
                _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);
                _mydata->auxDBFilename = "data/IsotropicLinearElasticityPlaneStrain_UniformStrain_aux.spatialdb";

            } // setUp

        }; // TestIsotropicLinearElasticityPlaneStrain_UniformStrain


        // ----------------------------------------------------------------------

        class TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP1 : public TestIsotropicLinearElasticityPlaneStrain_UniformStrain {

            CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP1,  TestIsotropicLinearElasticityPlaneStrain);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestIsotropicLinearElasticityPlaneStrain_UniformStrain::setUp();
                CPPUNIT_ASSERT(_mydata);

                _mydata->meshFilename = "data/tri_small.mesh";

                _mydata->numSolnFields = 1;
                static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
                    {1, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // disp
                };
                _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

                _initializeMin();
            } // setUp

        }; // TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP1
        CPPUNIT_TEST_SUITE_REGISTRATION(TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP1);


        // ----------------------------------------------------------------------
        class TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP2 : public TestIsotropicLinearElasticityPlaneStrain_UniformStrain {

            CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP2,  TestIsotropicLinearElasticityPlaneStrain);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestIsotropicLinearElasticityPlaneStrain_UniformStrain::setUp();
                CPPUNIT_ASSERT(_mydata);

                _mydata->meshFilename = "data/tri_small.mesh";

                _mydata->numSolnFields = 1;
                static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
                    {2, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // disp
                };
                _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

                static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
                    {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // density
                    {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // shear_modulus
                    {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // bulk_modulus
                };
                _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

                _initializeMin();
            } // setUp

        }; // TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP2
        CPPUNIT_TEST_SUITE_REGISTRATION(TestIsotropicLinearElasticityPlaneStrain_UniformStrain_TriP2);


        // ----------------------------------------------------------------------
        class TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadP1 : public TestIsotropicLinearElasticityPlaneStrain_UniformStrain { // TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadP1

            CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadP1,  TestIsotropicLinearElasticityPlaneStrain);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestIsotropicLinearElasticityPlaneStrain_UniformStrain::setUp();
                CPPUNIT_ASSERT(_mydata);

                _mydata->meshFilename = "data/quad_aligned.mesh";

                _mydata->numSolnFields = 1;
                static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
                    {1, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // disp
                };
                _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

                static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // density
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // shear_modulus
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // bulk_modulus
                };
                _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

                _initializeMin();
            } // setUp

        }; // TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadP1
        CPPUNIT_TEST_SUITE_REGISTRATION(TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadP1);


        // ----------------------------------------------------------------------
        class TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadP2 : public TestIsotropicLinearElasticityPlaneStrain_UniformStrain { // TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadP2

            CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadP2,  TestIsotropicLinearElasticityPlaneStrain);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestIsotropicLinearElasticityPlaneStrain_UniformStrain::setUp();
                CPPUNIT_ASSERT(_mydata);

                _mydata->meshFilename = "data/quad_aligned.mesh";

                _mydata->numSolnFields = 1;
                static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
                    {2, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // disp
                };
                _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

                static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
                    {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // density
                    {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // shear_modulus
                    {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // bulk_modulus
                };
                _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

                _initializeMin();
            } // setUp

        }; // TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadP2
        CPPUNIT_TEST_SUITE_REGISTRATION(TestIsotropicLinearElasticityPlaneStrain_UniformStrain_QuadP2);


    } // materials
} // pylith

// End of file
