// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestIsotropicLinearElasticityPlaneStrain.hh" // Implementation of cases

#include "pylith/materials/IsotropicLinearElasticityPlaneStrain.hh" // USES IsotropicLinearElasticityPlaneStrain
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// forward declarations
namespace pylith {
    namespace materials {
        class TestIsotropicLinearElasticityPlaneStrain_GravityRefState;

        class TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP1;
        class TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP2;
        class TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP3;
        class TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP4;

        class TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ1;
        class TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ2;
        class TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ3;
        class TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ4;

    } // materials
} // pylith

// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState :
    public pylith::materials::TestIsotropicLinearElasticityPlaneStrain {

    /// Spatial database user functions for auxiiliary subfields (includes derived fields).
    static const double SMALL;
    static const double GACC;
    static const double YMAX;

    // Density
    static double density(const double x,
                          const double y) {
        return 2500.0;
    } // density
    static const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    // Vs
    static double vs(const double x,
                     const double y) {
        return 3000.0;
    } // vs
    static const char* vs_units(void) {
        return "m/s";
    } // vs_units

    // Vp
    static double vp(const double x,
                     const double y) {
        return sqrt(3.0)*vs(x,y);
    } // vp
    static const char* vp_units(void) {
        return "m/s";
    } // vp_units

    // shear modulus
    static double shearModulus(const double x,
                               const double y) {
        return density(x,y) * vs(x,y) * vs(x,y);
    } // shearModulus
    static const char* shearModulus_units(void) {
        return "Pa";
    } // shearModulus_units

    // bulk modulus
    static double bulkModulus(const double x,
                              const double y) {
        return density(x,y)*(vp(x,y)*vp(x,y) - 4.0/3.0*vs(x,y)*vs(x,y));
    } // bulkModulus
    static const char* bulkModulus_units(void) {
        return "Pa";
    } // bulkModulus_units

    // Spatial database user functions for solution subfields.

    static double referenceMeanStress(const double x,
                                      const double y) {
        return -density(x,y) * GACC * (YMAX-y);
    } // referenceMeanStress
    static double referenceShearStress(const double x,
                                       const double y) {
        return 0.0;
    } // referenceShearStress
    static const char* stress_units(void) {
        return "Pa";
    } // stress_units
    static double referenceStrain(const double x,
                                  const double y) {
        return 0.0;
    } // referencStrain
    static const char* strain_units(void) {
        return "none";
    } // strain_units

    static double setGravityAcc_x(const double x,
                               const double y) {
        return 0.0;
    } // setGravityAcc_x
    static double setGravityAcc_y(const double x,
                               const double y) {
        return -GACC;
    } // setGravityAcc_y
    static const char* acc_units(void) {
        return "m/s**2";
    } // acc_units


    // Displacement
    static double disp_x(const double x,
                         const double y) {
        return 0.0;
    } // disp_x
    static double disp_y(const double x,
                         const double y) {
        return 0.0;
    } // disp_y
    static const char* disp_units(void) {
        return "m";
    } // disp_units

    static double disp_dot_x(const double x,
                             const double y) {
        return 0.0;
    } // disp_dot_x
    static double disp_dot_y(const double x,
                             const double y) {
        return 0.0;
    } // disp_dot_y
    static const char* disp_dot_units(void) {
        return "m/s";
    } // disp_dot_units


    // Displacement + perturbation
    static double disp_perturb_x(const double x,
                                 const double y) {
        return disp_x(x, y) + SMALL;
    } // disp_perturb_x
    static double disp_perturb_y(const double x,
                                 const double y) {
        return disp_y(x, y) + SMALL;
    } // disp_perturb_y

protected:

    void setUp(void) {
        TestIsotropicLinearElasticityPlaneStrain::setUp();
        _mydata = new TestIsotropicLinearElasticityPlaneStrain_Data(); CPPUNIT_ASSERT(_mydata);

        // dimension set in base class.
        // meshFilename set in derived class.
        _mydata->boundaryLabel = "boundary";

        CPPUNIT_ASSERT(_mydata->normalizer);
        _mydata->normalizer->setLengthScale(1.0e+03);
        _mydata->normalizer->setTimeScale(2.0);
        _mydata->normalizer->setPressureScale(2.25e+10);
        _mydata->normalizer->computeDensityScale();

        delete _mydata->gravityField; _mydata->gravityField = new spatialdata::spatialdb::GravityField();
        _mydata->gravityField->setGravityDir(0.0, -1.0, 0.0);
        _mydata->gravityField->setGravityAcc(GACC);

        _mydata->t = 1.0;
        _mydata->dt = 0.05;
        _mydata->s_tshift = 1.0 / _mydata->dt;

        // solnDiscretizations set in derived class.

        _mydata->numAuxSubfields = 6;
        static const char* _auxSubfields[6] = {"density", "shear_modulus", "bulk_modulus", "gravitational_acceleration", "reference_stress", "reference_strain" };
        _mydata->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
            pylith::topology::Field::Discretization(0, 1), // gravitational_acceleration
            pylith::topology::Field::Discretization(1, 1), // reference_stress
            pylith::topology::Field::Discretization(1, 1)  // reference_strain
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        CPPUNIT_ASSERT(_mydata->auxDB);
        _mydata->auxDB->addValue("density", density, density_units());
        _mydata->auxDB->addValue("vp", vp, vp_units());
        _mydata->auxDB->addValue("vs", vs, vs_units());
        _mydata->auxDB->addValue("shear_modulus", shearModulus, shearModulus_units());
        _mydata->auxDB->addValue("bulk_modulus", bulkModulus, bulkModulus_units());
        _mydata->auxDB->addValue("reference_stress_xx", referenceMeanStress, stress_units());
        _mydata->auxDB->addValue("reference_stress_yy", referenceMeanStress, stress_units());
        _mydata->auxDB->addValue("reference_stress_zz", referenceMeanStress, stress_units());
        _mydata->auxDB->addValue("reference_stress_xy", referenceShearStress, stress_units());
        _mydata->auxDB->addValue("reference_strain_xx", referenceStrain, strain_units());
        _mydata->auxDB->addValue("reference_strain_yy", referenceStrain, strain_units());
        _mydata->auxDB->addValue("reference_strain_zz", referenceStrain, strain_units());
        _mydata->auxDB->addValue("reference_strain_xy", referenceStrain, strain_units());
        _mydata->auxDB->addValue("gravitational_acceleration_x", setGravityAcc_x, acc_units()); // test of subfield.
        _mydata->auxDB->addValue("gravitational_acceleration_y", setGravityAcc_y, acc_units());

        CPPUNIT_ASSERT(_mydata->solnDB);
        _mydata->solnDB->addValue("displacement_x", disp_x, disp_units());
        _mydata->solnDB->addValue("displacement_y", disp_y, disp_units());
        _mydata->solnDB->addValue("displacement_dot_x", disp_dot_x, disp_dot_units());
        _mydata->solnDB->addValue("displacement_dot_y", disp_dot_y, disp_dot_units());

        CPPUNIT_ASSERT(_mydata->perturbDB);
        _mydata->perturbDB->addValue("displacement_x", disp_perturb_x, disp_units());
        _mydata->perturbDB->addValue("displacement_y", disp_perturb_y, disp_units());
        _mydata->perturbDB->addValue("displacement_dot_x", disp_dot_x, disp_dot_units());
        _mydata->perturbDB->addValue("displacement_dot_y", disp_dot_y, disp_dot_units());

        CPPUNIT_ASSERT(_mymaterial);
        _mymaterial->useInertia(false);
        _mymaterial->useBodyForce(false);
        _mymaterial->useReferenceState(true);

        _mymaterial->setLabel("Isotropic Linear Elascitity Plane Strain");
        _mymaterial->id(24);
    } // setUp

}; // TestIsotropicLinearElasticityPlaneStrain_GravityRefState
const double pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState::SMALL = 0.1;
const double pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState::GACC = 9.80665;
const double pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState::YMAX = +4.0e+3;

// ----------------------------------------------------------------------

class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP1 :
    public pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP1,
                           TestIsotropicLinearElasticityPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticityPlaneStrain_GravityRefState::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tri_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(1, 1)  // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP1);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP2 :
    public pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP2,
                           TestIsotropicLinearElasticityPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticityPlaneStrain_GravityRefState::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tri_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(2, 2), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
            pylith::topology::Field::Discretization(0, 2), // gravitational_acceleration
            pylith::topology::Field::Discretization(1, 2), // reference_stress
            pylith::topology::Field::Discretization(1, 2), // reference_strain
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP2);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP3 :
    public pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP3,
                           TestIsotropicLinearElasticityPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticityPlaneStrain_GravityRefState::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tri_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(3, 3), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
            pylith::topology::Field::Discretization(0, 3), // gravitational_acceleration
            pylith::topology::Field::Discretization(1, 3), // reference_stress
            pylith::topology::Field::Discretization(1, 3), // reference_strain
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP3);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP4 :
    public pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP4,
                           TestIsotropicLinearElasticityPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticityPlaneStrain_GravityRefState::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/tri_small.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(4, 4), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
            pylith::topology::Field::Discretization(0, 4), // gravitational_acceleration
            pylith::topology::Field::Discretization(1, 4), // reference_stress
            pylith::topology::Field::Discretization(1, 4), // reference_strain
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP4
// Leave this out for now to shorten runtime.
//CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState_TriP4);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ1 :
    public pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ1,
                           TestIsotropicLinearElasticityPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticityPlaneStrain_GravityRefState::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/quad_aligned.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(1, 1), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
            pylith::topology::Field::Discretization(0, 1), // gravitational_acceleration
            pylith::topology::Field::Discretization(1, 1), // reference_stress
            pylith::topology::Field::Discretization(1, 1), // reference_strain
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ1);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ2 :
    public pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ2,  TestIsotropicLinearElasticityPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticityPlaneStrain_GravityRefState::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/quad_aligned.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(2, 2), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
            pylith::topology::Field::Discretization(0, 2), // gravitational_acceleration
            pylith::topology::Field::Discretization(1, 2), // reference_stress
            pylith::topology::Field::Discretization(1, 2), // reference_strain
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ2);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ3 :
    public pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState {

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ3,
                           TestIsotropicLinearElasticityPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticityPlaneStrain_GravityRefState::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/quad_aligned.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(3, 3), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
            pylith::topology::Field::Discretization(0, 3), // gravitational_acceleration
            pylith::topology::Field::Discretization(1, 3), // reference_stress
            pylith::topology::Field::Discretization(1, 3), // reference_strain
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ3);


// ----------------------------------------------------------------------
class pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ4 : public pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState { // TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ4

    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ4,
                           TestIsotropicLinearElasticityPlaneStrain);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearElasticityPlaneStrain_GravityRefState::setUp();
        CPPUNIT_ASSERT(_mydata);

        _mydata->meshFilename = "data/quad_aligned.mesh";

        _mydata->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(4, 4), // disp
        };
        _mydata->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
            pylith::topology::Field::Discretization(0, 4), // gravitational_acceleration
            pylith::topology::Field::Discretization(1, 4), // reference_stress
            pylith::topology::Field::Discretization(1, 4), // reference_strain
        };
        _mydata->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _initializeMin();
    } // setUp

}; // TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ4
// Leave this out for now to shorten runtime.
//CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestIsotropicLinearElasticityPlaneStrain_GravityRefState_QuadQ4);


// End of file
