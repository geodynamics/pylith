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

#include "TestIsotropicLinearMaxwellPlaneStrain.hh" // Implementation of cases

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
namespace pylith {
    namespace materials {

        class TestIsotropicLinearMaxwellPlaneStrain_VarStrain : public TestIsotropicLinearMaxwellPlaneStrain {
protected:

            void setUp(void) {
                TestIsotropicLinearMaxwellPlaneStrain::setUp();
                _mydata = new TestIsotropicLinearMaxwellPlaneStrain_Data(); CPPUNIT_ASSERT(_mydata);

                // dimension set in base class.
                _mydata->useInertia = false;
                _mydata->useBodyForce = true;
                _mydata->useGravity = false;
                _mydata->useReferenceState = false;

                // meshFilename set in derived class.
                _mydata->materialLabel = "IsotropicLinearMaxwell";
                _mydata->materialId = 24;
                _mydata->boundaryLabel = "boundary";

                CPPUNIT_ASSERT(_mydata->normalizer);
                _mydata->normalizer->lengthScale(1.0e+03);
                _mydata->normalizer->timeScale(2.0);
                _mydata->normalizer->densityScale(3.0e+3);
                _mydata->normalizer->pressureScale(2.25e+10);

                _mydata->t = 1.0e7;
                _mydata->dt = 1.0e7;
                _mydata->tshift = 1.0 / _mydata->dt;

                // solnDiscretizations set in derived class.
                _mydata->solnDBFilename = "data/IsotropicLinearMaxwellPlaneStrain_VarStrain_soln.spatialdb";
                _mydata->pertDBFilename = "data/IsotropicLinearMaxwellPlaneStrain_VarStrain_pert.spatialdb";

                _mydata->numAuxFields = 7;
                static const char* _auxFields[7] = {"density", "shear_modulus", "bulk_modulus", "maxwell_time",
						    "total_strain", "viscous_strain", "body_force"};
                _mydata->auxFields = _auxFields;
                static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // density
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // shear_modulus
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // bulk_modulus
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // maxwell_time
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // total_strain
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // viscous_strain
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // body_force
                };
                _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);
                _mydata->auxDBFilename = "data/IsotropicLinearMaxwellPlaneStrain_VarStrain_aux.spatialdb";

            } // setUp

        }; // TestIsotropicLinearMaxwellPlaneStrain_VarStrain


        // ----------------------------------------------------------------------

        class TestIsotropicLinearMaxwellPlaneStrain_VarStrain_TriP1 : public TestIsotropicLinearMaxwellPlaneStrain_VarStrain {

            CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearMaxwellPlaneStrain_VarStrain_TriP1,  TestIsotropicLinearMaxwellPlaneStrain);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestIsotropicLinearMaxwellPlaneStrain_VarStrain::setUp();
                CPPUNIT_ASSERT(_mydata);

                _mydata->meshFilename = "data/tri_small.mesh";

                _mydata->numSolnFields = 1;
                static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
                    {1, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // disp
                };
                _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

                _initializeMin();
            } // setUp

        }; // TestIsotropicLinearMaxwellPlaneStrain_VarStrain_TriP1
        CPPUNIT_TEST_SUITE_REGISTRATION(TestIsotropicLinearMaxwellPlaneStrain_VarStrain_TriP1);


        // ----------------------------------------------------------------------
        class TestIsotropicLinearMaxwellPlaneStrain_VarStrain_TriP2 : public TestIsotropicLinearMaxwellPlaneStrain_VarStrain {

            CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearMaxwellPlaneStrain_VarStrain_TriP2,  TestIsotropicLinearMaxwellPlaneStrain);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestIsotropicLinearMaxwellPlaneStrain_VarStrain::setUp();
                CPPUNIT_ASSERT(_mydata);

                _mydata->meshFilename = "data/tri_small.mesh";

                _mydata->numSolnFields = 1;
                static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
                    {2, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // disp
                };
                _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

                static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
                    {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // density
                    {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // shear_modulus
                    {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // bulk_modulus
                    {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // maxwell_time
                    {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // total_strain
                    {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // viscous_strain
                    {0, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // body_force
                };
                _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

                _initializeMin();
            } // setUp

        }; // TestIsotropicLinearMaxwellPlaneStrain_VarStrain_TriP2
        CPPUNIT_TEST_SUITE_REGISTRATION(TestIsotropicLinearMaxwellPlaneStrain_VarStrain_TriP2);


        // ----------------------------------------------------------------------
        class TestIsotropicLinearMaxwellPlaneStrain_VarStrain_QuadP1 : public TestIsotropicLinearMaxwellPlaneStrain_VarStrain { // TestIsotropicLinearMaxwellPlaneStrain_VarStrain_QuadP1

            CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearMaxwellPlaneStrain_VarStrain_QuadP1,  TestIsotropicLinearMaxwellPlaneStrain);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestIsotropicLinearMaxwellPlaneStrain_VarStrain::setUp();
                CPPUNIT_ASSERT(_mydata);

                _mydata->meshFilename = "data/quad_aligned.mesh";

                _mydata->numSolnFields = 1;
                static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
                    {1, 2, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // disp
                };
                _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

                static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // density
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // shear_modulus
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // bulk_modulus
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // maxwell_time
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // total_strain
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // viscous_strain
                    {0, 1, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // body_force
                };
                _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

                _initializeMin();
            } // setUp

        }; // TestIsotropicLinearMaxwellPlaneStrain_VarStrain_QuadP1
        CPPUNIT_TEST_SUITE_REGISTRATION(TestIsotropicLinearMaxwellPlaneStrain_VarStrain_QuadP1);


        // ----------------------------------------------------------------------
        class TestIsotropicLinearMaxwellPlaneStrain_VarStrain_QuadP2 : public TestIsotropicLinearMaxwellPlaneStrain_VarStrain { // TestIsotropicLinearMaxwellPlaneStrain_VarStrain_QuadP2

            CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearMaxwellPlaneStrain_VarStrain_QuadP2,  TestIsotropicLinearMaxwellPlaneStrain);
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestIsotropicLinearMaxwellPlaneStrain_VarStrain::setUp();
                CPPUNIT_ASSERT(_mydata);

                _mydata->meshFilename = "data/quad_aligned.mesh";

                _mydata->numSolnFields = 1;
                static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
                    {2, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // disp
                };
                _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

                static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
                    {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // density
                    {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // shear_modulus
                    {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // bulk_modulus
                    {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // maxwell_time
                    {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // total_strain
                    {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // viscous_strain
                    {0, 3, true, pylith::topology::Field::POLYNOMIAL_SPACE}, // body_force
                };
                _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

                _initializeMin();
            } // setUp

        }; // TestIsotropicLinearMaxwellPlaneStrain_VarStrain_QuadP2
        CPPUNIT_TEST_SUITE_REGISTRATION(TestIsotropicLinearMaxwellPlaneStrain_VarStrain_QuadP2);


    } // materials
} // pylith

// End of file
